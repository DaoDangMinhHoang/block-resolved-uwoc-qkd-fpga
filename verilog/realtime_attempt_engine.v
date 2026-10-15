// ============================================================================
// realtime_attempt_engine.v   (rev F)
//
// Controlled-rate direct attempt engine for the Q2 security campaign.
//
// IMPORTANT:
// This path bypasses the already-validated OOK serializer/deserializer and
// exercises the stochastic uwoc_channel directly.
//
// At 50 MHz and rate_hz = 10_000_000:
//     one pacing request every 5 clocks.
//
// A minimum launch spacing of MIN_GAP clocks is enforced because uwoc_channel
// contains a pipelined h -> p_sig path. With MIN_GAP = 5, the next attempt is
// never launched before that pipeline has settled.
//
// The basis scheduler used here is deliberately pseudo-random rather than the
// original slow ring-oscillator TRNG. The original TRNG remains in the full
// OOK BB84 path. The realtime path is a high-rate channel-statistics mode.
//
// No missed request is queued or replayed. If a pacing deadline occurs while
// the engine cannot launch, overrun_count increments and that deadline is
// discarded. This prevents later bursts from hiding an unsustainable rate.
//
// Any paper run claiming real-time operation must satisfy:
//     timeout_count = 0
//     overrun_count = 0
//     attempts_per_sec ~= rate_hz
//
// [rev F, 2026-09-15]
//   * Basis schedulers are 64-bit xorshift (period 2^64-1) instead of 32-bit:
//     at 10 Mattempt/s a 32-bit state repeats after 2^32-1 attempts, i.e. every
//     65 536 blocks, shorter than the 286 000-block L1 capture.
//   * A third xorshift64 draws the intensity class of every attempt when
//     decoy_enable = 1: class 0 signal, 1 decoy, 2 vacuum, with probabilities
//     set by decoy_preset (table below). att_class is registered on the same
//     edge as sample_en and held for the whole attempt; evt_class travels with
//     evt_done / evt_lost. With decoy_enable = 0 every class is 0.
//   * All three generators are reseeded on `clear`, so the nominal-Y0 and
//     Y0 = 0 captures of a paired control see the same basis AND class
//     sequences.
// ============================================================================

module realtime_attempt_engine #(
    parameter integer CLK_HZ         = 50000000,
    parameter integer RESULT_TIMEOUT = 256,
    parameter integer MIN_GAP        = 5,

    parameter [63:0] A_BASIS_SEED =
        64'hA9C8_049F_12DA_AC6B,

    parameter [63:0] B_BASIS_SEED =
        64'hDACC_E005_EC05_AD8F,

    parameter [63:0] CLASS_SEED =
        64'h1B6F_01E1_67C4_752A
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire        run,
    input  wire        clear,

    // 0 = free run at the maximum safe rate.
    input  wire [31:0] rate_hz,

    // Registered response from uwoc_channel.
    input  wire        ch_click,
    input  wire        ch_no_click,
    input  wire        ch_err,

    // [rev F] Intensity classes for decoy-state campaigns.
    input  wire        decoy_enable,
    input  wire [1:0]  decoy_preset,

    // One-clock pulse per physical attempt.
    output reg         sample_en,

    // [rev F] Class of the attempt in flight, valid WITH sample_en.
    output reg  [1:0]  att_class,

    // One-clock event pulses.
    output reg         evt_done,
    output reg         evt_lost,
    output reg         evt_basis_match,
    output reg         evt_basis_is_z,
    output reg         evt_data_error,
    output reg  [1:0]  evt_class,

    output reg  [31:0] timeout_count,
    output reg  [31:0] overrun_count,

    output wire        busy
);

    // ------------------------------------------------------------------------
    // Parameter sizing for old Quartus releases
    // ------------------------------------------------------------------------

    localparam [15:0] TIMEOUT_U16 =
        RESULT_TIMEOUT;

    localparam [31:0] CLK_HZ_U32 =
        CLK_HZ;

    // If launches must be MIN_GAP clocks apart, a launch at clock n reloads
    // MIN_GAP-1. The counter reaches zero before clock n+MIN_GAP.
    //
    // Previous code loaded MIN_GAP itself, introducing an off-by-one:
    // MIN_GAP=5 rejected the pacing tick five clocks later and effectively
    // reduced a requested 10 MHz stream to about 5 MHz.
    localparam [7:0] GAP_RELOAD =
        (MIN_GAP > 0)
        ? (MIN_GAP - 1)
        : 0;


    localparam ST_IDLE =
        1'b0;

    localparam ST_WAIT =
        1'b1;


    reg state;

    reg [15:0] wait_cnt;

    reg [7:0] gap_cnt;


    // ------------------------------------------------------------------------
    // High-rate pseudo-random basis and class schedulers
    // ------------------------------------------------------------------------
    //
    // The original ring-oscillator TRNG cannot provide a fresh independent
    // basis bit at 10 Mattempt/s. It therefore remains on the full OOK path.
    //
    // For the realtime measurement path we use deterministic high-rate PRNGs
    // (64-bit xorshift, Marsaglia 2003). Each state advances exactly once per
    // ACTUAL launched attempt, and the most significant bit is used.
    // ------------------------------------------------------------------------

    reg [63:0] prng_a;
    reg [63:0] prng_b;
    reg [63:0] prng_c;

    reg a_basis_q;
    reg b_basis_q;


    function [63:0] xorshift64;

        input [63:0] state_in;

        reg [63:0] x;

        begin

            x = state_in;

            x = x ^ (x << 13);

            x = x ^ (x >> 7);

            x = x ^ (x << 17);

            xorshift64 = x;

        end
    endfunction


    // Class probabilities, 16-bit thresholds on prng_c[63:48]:
    //   preset 0: signal 1/2, decoy 1/4, vacuum 1/4
    //   preset 1: signal 0.70, decoy 0.20, vacuum 0.10
    //   preset 2: signal 0.80, decoy 0.15, vacuum 0.05
    //   preset 3: signal 1/3, decoy 1/3, vacuum 1/3
    reg [15:0] th_signal;
    reg [15:0] th_decoy;

    always @(*) begin

        case (decoy_preset)

            2'd0: begin
                th_signal = 16'd32768;
                th_decoy  = 16'd49152;
            end

            2'd1: begin
                th_signal = 16'd45875;
                th_decoy  = 16'd58982;
            end

            2'd2: begin
                th_signal = 16'd52429;
                th_decoy  = 16'd62259;
            end

            default: begin
                th_signal = 16'd21845;
                th_decoy  = 16'd43690;
            end

        endcase
    end


    wire [15:0] class_word =
        prng_c[63:48];

    wire [1:0] class_next =
        (!decoy_enable)            ? 2'd0 :
        (class_word < th_signal)   ? 2'd0 :
        (class_word < th_decoy)    ? 2'd1 :
                                     2'd2;


    wire bmatch_q =
        (
            a_basis_q ==
            b_basis_q
        );


    // Alice encoding:
    // 0 = Z
    // 1 = X
    wire basis_is_z_q =
        ~a_basis_q;


    // ------------------------------------------------------------------------
    // Rate pacer
    // ------------------------------------------------------------------------

    reg [32:0] acc;

    reg tick;


    wire free_run =
        (
            rate_hz ==
            32'd0
        );


    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin

            acc <=
                33'd0;

            tick <=
                1'b0;

        end

        else if (
            clear ||
            !run ||
            free_run
        ) begin

            acc <=
                33'd0;

            tick <=
                1'b0;

        end

        else begin

            if (
                acc +
                {1'b0, rate_hz}
                >=
                {1'b0, CLK_HZ_U32}
            ) begin

                acc <=
                    acc +
                    {1'b0, rate_hz}
                    -
                    {1'b0, CLK_HZ_U32};

                tick <=
                    1'b1;

            end

            else begin

                acc <=
                    acc +
                    {1'b0, rate_hz};

                tick <=
                    1'b0;

            end
        end
    end


    // ------------------------------------------------------------------------
    // Attempt sequencer
    // ------------------------------------------------------------------------

    assign busy =
        (
            state !=
            ST_IDLE
        );


    wire result_valid =
        ch_click |
        ch_no_click;


    wire can_launch =
        (
            state ==
            ST_IDLE
        )
        &&
        (
            gap_cnt ==
            8'd0
        );


    always @(posedge clk or negedge rst_n) begin

        if (!rst_n) begin

            state <=
                ST_IDLE;

            wait_cnt <=
                16'd0;

            gap_cnt <=
                8'd0;


            prng_a <=
                A_BASIS_SEED;

            prng_b <=
                B_BASIS_SEED;

            prng_c <=
                CLASS_SEED;


            a_basis_q <=
                1'b0;

            b_basis_q <=
                1'b0;

            att_class <=
                2'd0;


            sample_en <=
                1'b0;


            evt_done <=
                1'b0;

            evt_lost <=
                1'b0;

            evt_basis_match <=
                1'b0;

            evt_basis_is_z <=
                1'b0;

            evt_data_error <=
                1'b0;

            evt_class <=
                2'd0;


            timeout_count <=
                32'd0;

            overrun_count <=
                32'd0;

        end

        else if (clear) begin

            state <=
                ST_IDLE;

            wait_cnt <=
                16'd0;

            gap_cnt <=
                8'd0;


            // Same basis and class sequences are regenerated for each
            // campaign. This is useful when comparing nominal Y0 against Y0=0.
            prng_a <=
                A_BASIS_SEED;

            prng_b <=
                B_BASIS_SEED;

            prng_c <=
                CLASS_SEED;


            a_basis_q <=
                1'b0;

            b_basis_q <=
                1'b0;

            att_class <=
                2'd0;


            sample_en <=
                1'b0;


            evt_done <=
                1'b0;

            evt_lost <=
                1'b0;

            evt_basis_match <=
                1'b0;

            evt_basis_is_z <=
                1'b0;

            evt_data_error <=
                1'b0;

            evt_class <=
                2'd0;


            timeout_count <=
                32'd0;

            overrun_count <=
                32'd0;

        end

        else begin

            // Default all pulse outputs low.
            sample_en <=
                1'b0;

            evt_done <=
                1'b0;

            evt_lost <=
                1'b0;

            evt_basis_match <=
                1'b0;

            evt_basis_is_z <=
                1'b0;

            evt_data_error <=
                1'b0;


            // Minimum inter-attempt spacing.
            if (
                gap_cnt !=
                8'd0
            ) begin

                gap_cnt <=
                    gap_cnt -
                    1'b1;

            end


            if (!run) begin

                state <=
                    ST_IDLE;

                wait_cnt <=
                    16'd0;

                gap_cnt <=
                    8'd0;

            end

            else begin

                // ------------------------------------------------------------
                // Deadline integrity check
                // ------------------------------------------------------------
                //
                // A paced request arriving while either:
                //   - the previous attempt is still in flight, OR
                //   - the minimum pipeline spacing has not expired
                //
                // is a missed real-time deadline.
                //
                // The old implementation counted only ST_WAIT and silently lost
                // ticks that arrived while ST_IDLE but gap_cnt != 0.
                // ------------------------------------------------------------

                if (
                    !free_run &&
                    tick &&
                    !can_launch
                ) begin

                    overrun_count <=
                        overrun_count +
                        1'b1;

                end


                case (state)

                    // ========================================================
                    ST_IDLE: begin

                        if (
                            can_launch &&
                            (
                                free_run ||
                                tick
                            )
                        ) begin

                            // Freeze a deterministic pseudo-random basis pair
                            // and intensity class for this complete attempt.
                            a_basis_q <=
                                prng_a[63];

                            b_basis_q <=
                                prng_b[63];

                            att_class <=
                                class_next;


                            // Advance exactly once per ACTUAL attempt.
                            prng_a <=
                                xorshift64(
                                    prng_a
                                );

                            prng_b <=
                                xorshift64(
                                    prng_b
                                );

                            prng_c <=
                                xorshift64(
                                    prng_c
                                );


                            sample_en <=
                                1'b1;

                            wait_cnt <=
                                16'd0;

                            gap_cnt <=
                                GAP_RELOAD;

                            state <=
                                ST_WAIT;

                        end
                    end


                    // ========================================================
                    ST_WAIT: begin

                        if (result_valid) begin

                            evt_done <=
                                ch_click;

                            evt_lost <=
                                ch_no_click;

                            // Basis fields remain valid even on no-click.
                            evt_basis_match <=
                                bmatch_q;

                            evt_basis_is_z <=
                                basis_is_z_q;

                            evt_data_error <=
                                bmatch_q &&
                                ch_click &&
                                ch_err;

                            evt_class <=
                                att_class;

                            state <=
                                ST_IDLE;

                        end

                        else if (
                            wait_cnt >=
                            TIMEOUT_U16
                        ) begin

                            // Do NOT convert a hardware timeout into a plausible
                            // physical no-click. No attempt event is emitted.
                            //
                            // uwoc_channel has already advanced qub_index, so the
                            // next valid event creates index_step_bad in the
                            // security monitor and invalidates the capture.
                            timeout_count <=
                                timeout_count +
                                1'b1;

                            state <=
                                ST_IDLE;

                        end

                        else begin

                            wait_cnt <=
                                wait_cnt +
                                1'b1;

                        end
                    end


                    default: begin

                        state <=
                            ST_IDLE;

                    end

                endcase
            end
        end
    end

endmodule
