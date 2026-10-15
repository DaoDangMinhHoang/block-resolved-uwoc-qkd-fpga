// ============================================================
// CHANNEL MONITOR v3 — for the UWOC channel (replaces v2, used for atmospheric FSO)
// ============================================================
//
// WHY THE REWRITE? (quantitative conclusions from python/uwoc_channel_model.py)
//
//  [1] A 256-ATTEMPT WINDOW IS FAR TOO SMALL.
//      Underwater P_click ~ 10⁻³…10⁻⁴ (atmosphere: ~O(1)). A 256-attempt window
//      collects only 0.4 clicks → the window QBER is pure SHOT NOISE, not
//      information about the channel. CHECK 7 shows std(QBER) at WEAK
//      turbulence (24.4%) is LARGER than at SEVERE (22.7%) — i.e. the
//      controller would react backwards. → window raised to 2^ATTEMPT_LOG2
//      (default 16384) plus a new window_valid flag.
//
//  [2] THE OLD SNR PROXY MIXED STATIC LOSS WITH DYNAMIC FADING.
//      v2 derived snr from (qubits received) + (duty ratio of signal_detect).
//      Underwater, qubits received are dominated by exp(−c·d) — a constant
//      of the link configuration, NOT something the controller can improve.
//      → normalise by expectation: snr_level = 128·(actual / expected clicks),
//        using nexp_inv from the uwoc_channel ROM (multiply instead of divide).
//      Drop the signal_detect input: OOK duty cycle carries no channel info.
//
//  [3] NO TURBULENCE INDICATOR.
//      CHECK 8: across 5 turbulence levels the MEAN QBER moves only 70% while
//      the STANDARD DEVIATION between windows moves 252%. Underwater turbulence
//      shows up in VARIABILITY, not in the mean (P_click ≈ linear in h, so the
//      expectation cancels the fading). → add qber_jitter (EWMA of |ΔQBER|).
//
//  [4] SEPARATE OUT PHOTON LOSS. The dominant underwater event is no-click,
//      not a bit error → add loss_rate.
//
// OUTPUT SCALING
//   qber        0–200, one unit = 0.5%      (unchanged from v2)
//   snr_level   0–255, 128 = EXACTLY the link's nominal click rate
//   loss_rate   0–255, 255 = everything lost
//   qber_jitter 0–255, same scale as qber (0.5%/unit), an EWMA of |ΔQBER|
// ============================================================

module channel_monitor #(
    // Monitoring window = 2^ATTEMPT_LOG2 attempts.
    //
    // TRADE-OFF (CHECK 7 of uwoc_channel_model.py): the window is SQUEEZED from
    // both sides — too short and QBER drowns in shot noise, too long and fading
    // is averaged away so the controller no longer sees turbulence. With the
    // default link parameters, 2^16 = 65536 attempts give ~40–170 clicks at the
    // working range — enough samples while still inside one coherence block.
    parameter ATTEMPT_LOG2 = 16,

    // log2 of the --window parameter used to generate nexp_inv (uwoc_lut_gen.py).
    // Kept separate from ATTEMPT_LOG2 so simulation can shorten the window WITHOUT
    // regenerating the ROM; the difference is compensated by snr_shift.
    parameter NEXP_LOG2 = 16,

    // Minimum sifted qubits before the QBER estimate is considered trustworthy.
    // 16 samples → standard error of QBER ≈ 5% absolute at QBER 5%: enough to
    // separate the 4/8/11/15% thresholds. Setting it higher (32) makes a 15 m
    // clear-ocean link — perfectly healthy — get flagged as invalid.
    parameter MIN_SIFT = 16
) (
    input  wire        clk,
    input  wire        rst_n,

    // Qubit-level events (from the top_module FSM)
    input  wire        evt_qubit_done,   // click occurred
    input  wire        evt_qubit_lost,   // NO click (photon lost)
    input  wire        evt_basis_match,
    input  wire        evt_data_error,

    // Expected clicks, from the uwoc_channel ROM: 32768 / E[clicks per window]
    input  wire [15:0] nexp_inv,

    input  wire        enable,
    input  wire        clear,

    output reg  [7:0]  qber,          // 0–200 (0.5%/unit)
    output reg  [7:0]  snr_level,     // NORMALISED, 128 = nominal
    output reg  [7:0]  photon_rate,   // clicks/window, SATURATES at 255
    output reg  [15:0] photon_count,  // clicks/window, does NOT saturate.
                                      // Needed by the wavelength hill-climber:
                                      // it compares click rate between two λ,
                                      // and at short range photon_rate is
                                      // clamped at 255, killing the comparison.
    output reg  [7:0]  sifted_rate,
    output reg  [7:0]  loss_rate,     // 255 = everything lost
    output reg  [7:0]  qber_jitter,   // EWMA |ΔQBER| — turbulence indicator
    output reg         window_valid,  // enough samples to trust the estimate
    output reg         window_pulse
);

    localparam integer AW = ATTEMPT_LOG2;

    reg [AW-1:0] w_attempt;
    reg [15:0]   w_received, w_sifted, w_errors, w_lost;
    reg [7:0]    qber_prev;

    wire window_end = enable & (evt_qubit_done | evt_qubit_lost)
                    & (&w_attempt);

    // ---- Normalised SNR: 128 · (actual clicks / expected clicks) ----
    // nexp_inv = 32768 / E[clicks in a 2^NEXP_LOG2 window]
    //   ⇒ received·nexp_inv = 32768 at exactly nominal ⇒ >>8 gives 128.
    // If the real window is shorter than the LUT window, received shrinks by
    // exactly 2^(NEXP_LOG2−ATTEMPT_LOG2) → compensate by shifting right LESS.
    // One 16×16 multiply, NO divider needed.
    localparam integer SNR_SHIFT = 8 - (NEXP_LOG2 - ATTEMPT_LOG2);

    wire [31:0] snr_prod  = w_received * nexp_inv;
    wire [31:0] snr_shift = snr_prod >> SNR_SHIFT;
    wire [7:0]  snr_calc  = (snr_shift > 32'd255) ? 8'd255 : snr_shift[7:0];

    // ---- loss_rate: w_lost mapped onto the window's 0–255 scale ----
    wire [15:0] loss_sh = (AW >= 8) ? (w_lost >> (AW - 8)) : (w_lost << (8 - AW));
    wire [7:0]  loss_calc = (loss_sh > 16'd255) ? 8'd255 : loss_sh[7:0];

    wire [7:0] photon_calc = (w_received > 16'd255) ? 8'd255 : w_received[7:0];
    wire [7:0] sifted_calc = (w_sifted   > 16'd255) ? 8'd255 : w_sifted[7:0];

    wire [7:0] qber_calc = (w_sifted == 16'd0) ? 8'd200
                                               : approx_qber(w_errors, w_sifted);

    // ---- jitter: EWMA with coefficient 1/4 of |QBER(k) − QBER(k−1)| ----
    wire [7:0] qdiff = (qber_calc > qber_prev) ? (qber_calc - qber_prev)
                                               : (qber_prev - qber_calc);
    wire [7:0] jit_next = qber_jitter - (qber_jitter >> 2) + (qdiff >> 2);

    // ⚠ `clear` MUST sit in its own branch, it may NOT be merged into
    //   `if (!rst_n || clear)`. In an asynchronous reset branch Quartus only
    //   accepts signals that APPEAR in the always sensitivity list:
    //     Error (10200): cannot match operand(s) in the condition to the
    //                    corresponding edges in the enclosing event control
    //   Merged, the whole block is treated as combinational → latches inferred
    //   for w_attempt, qber, snr_level… then Error (12152): ch_mon fails to elaborate.
    //   ModelSim compiles it happily, so the failure only surfaces in quartus_map.
    //   Splitting as below preserves behaviour: `clear` was already SYNCHRONOUS
    //   (absent from the sensitivity list, it only acts on the clk rising edge).
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            w_attempt    <= {AW{1'b0}};
            w_received   <= 16'd0;
            w_sifted     <= 16'd0;
            w_errors     <= 16'd0;
            w_lost       <= 16'd0;
            qber         <= 8'd0;
            qber_prev    <= 8'd0;
            qber_jitter  <= 8'd0;
            snr_level    <= 8'd128;
            photon_rate  <= 8'd0;
            photon_count <= 16'd0;
            sifted_rate  <= 8'd0;
            loss_rate    <= 8'd0;
            window_valid <= 1'b0;
            window_pulse <= 1'b0;
        end else if (clear) begin
            w_attempt    <= {AW{1'b0}};
            w_received   <= 16'd0;
            w_sifted     <= 16'd0;
            w_errors     <= 16'd0;
            w_lost       <= 16'd0;
            qber         <= 8'd0;
            qber_prev    <= 8'd0;
            qber_jitter  <= 8'd0;
            snr_level    <= 8'd128;
            photon_rate  <= 8'd0;
            photon_count <= 16'd0;
            sifted_rate  <= 8'd0;
            loss_rate    <= 8'd0;
            window_valid <= 1'b0;
            window_pulse <= 1'b0;
        end else begin
            window_pulse <= 1'b0;

            if (enable) begin
                if (window_end) begin
                    // ---- Latch the window results ----
                    window_pulse <= 1'b1;
                    qber         <= qber_calc;
                    qber_prev    <= qber_calc;
                    qber_jitter  <= jit_next;
                    snr_level    <= snr_calc;
                    photon_rate  <= photon_calc;
                    photon_count <= w_received;
                    sifted_rate  <= sifted_calc;
                    loss_rate    <= loss_calc;
                    window_valid <= (w_sifted >= MIN_SIFT);

                    // ---- Open a new window, counting the current event ----
                    w_attempt  <= {AW{1'b0}};
                    w_received <= evt_qubit_done ? 16'd1 : 16'd0;
                    w_lost     <= evt_qubit_lost ? 16'd1 : 16'd0;
                    w_sifted   <= (evt_qubit_done & evt_basis_match) ? 16'd1 : 16'd0;
                    w_errors   <= (evt_qubit_done & evt_basis_match
                                                  & evt_data_error) ? 16'd1 : 16'd0;
                end else begin
                    if (evt_qubit_done) begin
                        w_attempt  <= w_attempt  + 1'b1;
                        w_received <= w_received + 1'b1;
                        if (evt_basis_match) begin
                            w_sifted <= w_sifted + 1'b1;
                            if (evt_data_error)
                                w_errors <= w_errors + 1'b1;
                        end
                    end
                    if (evt_qubit_lost) begin
                        w_attempt <= w_attempt + 1'b1;
                        w_lost    <= w_lost    + 1'b1;
                    end
                end
            end
        end
    end

    // ============================================================
    // QBER approximated by compare-and-shift (no divider)
    // 12 levels: 0, 1, 2, 4, 6, 8, 12, 16, 25, 33, 50, 100 (%)
    // Widened from 8-bit (v2) to 16-bit: underwater w_sifted can reach several
    // thousand within a 16384-attempt window.
    // ============================================================
    function [7:0] approx_qber;
        input [15:0] errors, sifted;
        reg [23:0] e1, e3, e6, e12;
    begin
        // Defaults so every branch drives the temporaries. Only the else
        // branch below actually uses them, but Quartus treats function
        // variables as static and warns about latch behaviour (10776)
        // when a path leaves them unassigned.
        e1  = 24'd0;
        e3  = 24'd0;
        e6  = 24'd0;
        e12 = 24'd0;

        if (errors == 16'd0)
            approx_qber = 8'd0;
        else if (errors >= sifted)
            approx_qber = 8'd200;                    // ≥ 100%
        else begin
            e1  = {8'd0, errors};
            e3  = e1 + e1 + e1;
            e6  = e3 + e3;
            e12 = e6 + e6;
            if      ({errors, 1'b0}     >= {8'd0, sifted}) approx_qber = 8'd100; // 50%
            else if (e3                 >= {8'd0, sifted}) approx_qber = 8'd66;  // 33%
            else if ({errors, 2'b00}    >= {8'd0, sifted}) approx_qber = 8'd50;  // 25%
            else if (e6                 >= {8'd0, sifted}) approx_qber = 8'd32;  // 16%
            else if ({errors, 3'b000}   >= {8'd0, sifted}) approx_qber = 8'd24;  // 12%
            else if (e12                >= {8'd0, sifted}) approx_qber = 8'd16;  // 8%
            else if ({errors, 4'd0}     >= {8'd0, sifted}) approx_qber = 8'd12;  // 6%
            else if ({errors, 5'd0}     >= {8'd0, sifted}) approx_qber = 8'd8;   // 4%
            else if ({errors, 6'd0}     >= {8'd0, sifted}) approx_qber = 8'd4;   // 2%
            else if ({errors, 7'd0}     >= {8'd0, sifted}) approx_qber = 8'd2;   // 1%
            else                                           approx_qber = 8'd1;   // <1%
        end
    end
    endfunction

endmodule
