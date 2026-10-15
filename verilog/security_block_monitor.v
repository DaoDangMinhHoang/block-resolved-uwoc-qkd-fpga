// ============================================================================
// security_block_monitor.v   (rev E)
//
// Block-resolved affine statistics, basis-resolved, for the journal extension.
// Runs in parallel with channel_monitor; the adaptive controller is untouched.
//
// Per block of exactly 2^ATTEMPT_LOG2 attempts:
//   attempts                      physical rounds  -> the mixture weight p_j
//   trials_z/x                    basis-matched attempts (click or lost)
//   clicks_z/x, errors_z/x        basis-matched detections and errors
//   clicks_all, lost_all          every attempt, matched or not
//   h_mean/h_min/h_max            fading of the block, Q8 (256 = 1.0)
//
// [rev E, 2026-09-15] Intensity classes for decoy-state campaigns.
//   The totals above still count EVERY attempt, whatever its class, so a run
//   with decoys disabled produces exactly the rev D numbers. On top of them the
//   monitor now keeps the same eight counters for class 1 (decoy) and class 2
//   (vacuum), packed 18 bits per field:
//       [ 17:  0] attempts   [ 35: 18] clicks_all
//       [ 53: 36] trials_z   [ 71: 54] clicks_z   [ 89: 72] errors_z
//       [107: 90] trials_x   [125:108] clicks_x   [143:126] errors_x
//   Class 0 (signal) is the total minus the two, computed on the host.
//   cfg_bus grows to 80 bits; bits [79:72] carry the decoy configuration so a
//   decoy change inside a block raises config_changed like any other change.
//
// Offline: G = clicks/attempts, T = errors/attempts, weight = attempts.
// The mixture weight in Proposition 1 is the number of physical rounds, so the
// denominator and the weight must both be `attempts`. Using trials as the
// weight folds basis-sampling noise into p_j for no reason.
// QBER is never formed here: it is a ratio, not an affine statistic.
//
// ---------------------------------------------------------------------------
// DESIGN INVARIANT -- read before changing the attempt engine.
//
// The window boundary comes from uwoc_channel's own attempt counter:
// `qub_index <= qub_cnt` latches the index of the qubit being sampled at that
// edge, so
//     block(this qubit) == qub_index >> (coh_sel + 5)
// holds with no off-by-one against the h that was applied. Setting
// coh_sel = ATTEMPT_LOG2 - 5 therefore puts exactly one fading realisation in
// every security block, and h_min == h_max becomes an assertion, not a filter.
//
// This holds ONLY while the attempt generator is strictly serial: it must not
// pulse sample_en for attempt n+1 before the events of attempt n have been
// emitted, or qub_index will already have advanced. realtime_attempt_engine.v
// is serial by construction. index_step_bad below is the runtime guard.
// ---------------------------------------------------------------------------
//
// cfg_bus packing (identical in top_module and in the Python decoder):
//   [ 1: 0] water       [ 5: 2] dist_idx   [ 8: 6] turb_level
//   [10: 9] lambda_sel  [14:11] mu_level   [16:15] adapt_mode
//   [40:17] p_noise     [44:41] coh_sel    [68:45] slot_width
//   [71:69] rate_sel    -- so a pacing change mid-block also sets config_changed
//   [72]    decoy_enable  [74:73] decoy_shift code  [76:75] decoy preset
//   [78:77] mu_gain_log2  [79] reserved
// ============================================================================

module security_block_monitor #(
    parameter ATTEMPT_LOG2 = 16
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire        enable,
    input  wire        clear,

    input  wire        evt_qubit_done,
    input  wire        evt_qubit_lost,
    input  wire        evt_basis_match,     // valid on lost attempts too
    input  wire        evt_basis_is_z,      // valid on lost attempts too
    input  wire        evt_data_error,
    input  wire [1:0]  evt_class,           // [rev E] 0 signal, 1 decoy, 2 vacuum

    input  wire [11:0] h_f_sample,
    input  wire [23:0] qub_index,
    input  wire [79:0] cfg_bus,

    output reg         block_pulse,
    output reg         block_synced,        // valid WITH block_pulse
    output reg  [31:0] block_id,

    output reg  [17:0] attempts,
    output reg  [17:0] clicks_all,
    output reg  [17:0] lost_all,
    output reg  [17:0] trials_z,
    output reg  [17:0] clicks_z,
    output reg  [17:0] errors_z,
    output reg  [17:0] trials_x,
    output reg  [17:0] clicks_x,
    output reg  [17:0] errors_x,

    output reg  [143:0] decoy_counts,       // [rev E] class 1, packed as above
    output reg  [143:0] vacuum_counts,      // [rev E] class 2, packed as above

    output reg  [11:0] h_mean,
    output reg  [11:0] h_min,
    output reg  [11:0] h_max,

    output reg  [79:0] blk_cfg,
    output reg         config_changed,
    output reg  [31:0] protocol_errors
);

    wire evt_attempt = evt_qubit_done | evt_qubit_lost;
    wire bad_attempt = evt_qubit_done & evt_qubit_lost;

    wire [ATTEMPT_LOG2-1:0] w_pos = qub_index[ATTEMPT_LOG2-1:0];
    wire window_last  = &w_pos;
    wire window_first = ~|w_pos;

    wire z_hit = evt_basis_match &  evt_basis_is_z;
    wire x_hit = evt_basis_match & ~evt_basis_is_z;

    wire [17:0] i_att   = 18'd1;
    wire [17:0] i_click = evt_qubit_done ? 18'd1 : 18'd0;
    wire [17:0] i_lost  = evt_qubit_lost ? 18'd1 : 18'd0;
    wire [17:0] i_tz    = z_hit ? 18'd1 : 18'd0;
    wire [17:0] i_cz    = (z_hit & evt_qubit_done) ? 18'd1 : 18'd0;
    wire [17:0] i_ez    = (z_hit & evt_qubit_done & evt_data_error) ? 18'd1 : 18'd0;
    wire [17:0] i_tx    = x_hit ? 18'd1 : 18'd0;
    wire [17:0] i_cx    = (x_hit & evt_qubit_done) ? 18'd1 : 18'd0;
    wire [17:0] i_ex    = (x_hit & evt_qubit_done & evt_data_error) ? 18'd1 : 18'd0;

    // [rev E] Increment vector of one class, same eight fields as the totals.
    function [143:0] class_inc;
        input hit;
        input done;
        input zh;
        input xh;
        input derr;
        begin
            class_inc = 144'd0;
            if (hit) begin
                class_inc[  0 +: 18] = 18'd1;
                class_inc[ 18 +: 18] = done               ? 18'd1 : 18'd0;
                class_inc[ 36 +: 18] = zh                 ? 18'd1 : 18'd0;
                class_inc[ 54 +: 18] = (zh & done)        ? 18'd1 : 18'd0;
                class_inc[ 72 +: 18] = (zh & done & derr) ? 18'd1 : 18'd0;
                class_inc[ 90 +: 18] = xh                 ? 18'd1 : 18'd0;
                class_inc[108 +: 18] = (xh & done)        ? 18'd1 : 18'd0;
                class_inc[126 +: 18] = (xh & done & derr) ? 18'd1 : 18'd0;
            end
        end
    endfunction

    // Field-wise addition of two packed class vectors.
    function [143:0] class_add;
        input [143:0] a;
        input [143:0] b;
        integer k;
        begin
            class_add = 144'd0;
            for (k = 0; k < 8; k = k + 1)
                class_add[k*18 +: 18] = a[k*18 +: 18] + b[k*18 +: 18];
        end
    endfunction

    wire [143:0] i_dec = class_inc(evt_class == 2'd1, evt_qubit_done,
                                   z_hit, x_hit, evt_data_error);
    wire [143:0] i_vac = class_inc(evt_class == 2'd2, evt_qubit_done,
                                   z_hit, x_hit, evt_data_error);

    reg [17:0] c_att, c_click, c_lost;
    reg [17:0] c_tz, c_cz, c_ez, c_tx, c_cx, c_ex;
    reg [143:0] c_dec, c_vac;
    reg [27:0] c_hsum;
    reg [11:0] c_hmin, c_hmax;
    reg [79:0] c_cfg;
    reg        c_cfg_chg;
    reg        armed;
    reg [31:0] next_block_id;
    reg [23:0] prev_index;
    reg        have_prev;

    // h_min/h_max are seeded at 0xFFF / 0x000 so the first compare of a fresh
    // block adopts the sample with no special priming branch.
    wire [27:0] hsum_n = c_hsum + {16'd0, h_f_sample};
    wire [11:0] hmin_n = (h_f_sample < c_hmin) ? h_f_sample : c_hmin;
    wire [11:0] hmax_n = (h_f_sample > c_hmax) ? h_f_sample : c_hmax;

    wire cfg_mismatch   = (cfg_bus != c_cfg);
    wire index_step_bad = have_prev &&
                          (qub_index != ((prev_index + 24'd1) & 24'hFFFFFF));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            block_pulse     <= 1'b0;
            block_synced    <= 1'b0;
            block_id        <= 32'd0;
            attempts        <= 18'd0;
            clicks_all      <= 18'd0;
            lost_all        <= 18'd0;
            trials_z        <= 18'd0;
            clicks_z        <= 18'd0;
            errors_z        <= 18'd0;
            trials_x        <= 18'd0;
            clicks_x        <= 18'd0;
            errors_x        <= 18'd0;
            decoy_counts    <= 144'd0;
            vacuum_counts   <= 144'd0;
            h_mean          <= 12'd0;
            h_min           <= 12'd0;
            h_max           <= 12'd0;
            blk_cfg         <= 80'd0;
            config_changed  <= 1'b0;
            protocol_errors <= 32'd0;
            armed           <= 1'b0;
            next_block_id   <= 32'd0;
            prev_index      <= 24'd0;
            have_prev       <= 1'b0;
            c_cfg           <= 80'd0;
            c_att   <= 18'd0; c_click <= 18'd0; c_lost <= 18'd0;
            c_tz    <= 18'd0; c_cz    <= 18'd0; c_ez   <= 18'd0;
            c_tx    <= 18'd0; c_cx    <= 18'd0; c_ex   <= 18'd0;
            c_dec   <= 144'd0; c_vac  <= 144'd0;
            c_hsum  <= 28'd0; c_hmin  <= 12'hFFF; c_hmax <= 12'd0;
            c_cfg_chg <= 1'b0;
        end
        else if (clear) begin
            // uwoc_channel restarts qub_index on the same pulse (stats_rst), so
            // both grids realign here.
            block_pulse     <= 1'b0;
            block_synced    <= 1'b0;
            protocol_errors <= 32'd0;
            armed           <= 1'b0;
            next_block_id   <= 32'd0;
            have_prev       <= 1'b0;
            c_att   <= 18'd0; c_click <= 18'd0; c_lost <= 18'd0;
            c_tz    <= 18'd0; c_cz    <= 18'd0; c_ez   <= 18'd0;
            c_tx    <= 18'd0; c_cx    <= 18'd0; c_ex   <= 18'd0;
            c_dec   <= 144'd0; c_vac  <= 144'd0;
            c_hsum  <= 28'd0; c_hmin  <= 12'hFFF; c_hmax <= 12'd0;
            c_cfg_chg <= 1'b0;
        end
        else begin
            block_pulse  <= 1'b0;
            block_synced <= 1'b0;

            if (enable && evt_attempt) begin
                prev_index <= qub_index;
                have_prev  <= 1'b1;

                if (bad_attempt || index_step_bad)
                    protocol_errors <= protocol_errors + 32'd1;

                if (!armed) begin
                    // Never report a partial first block: wait for a real
                    // window boundary before accumulating anything.
                    if (window_first) begin
                        armed     <= 1'b1;
                        c_att     <= i_att;   c_click <= i_click; c_lost <= i_lost;
                        c_tz      <= i_tz;    c_cz    <= i_cz;    c_ez   <= i_ez;
                        c_tx      <= i_tx;    c_cx    <= i_cx;    c_ex   <= i_ex;
                        c_dec     <= i_dec;   c_vac   <= i_vac;
                        c_hsum    <= {16'd0, h_f_sample};
                        c_hmin    <= h_f_sample;
                        c_hmax    <= h_f_sample;
                        c_cfg     <= cfg_bus;
                        c_cfg_chg <= 1'b0;
                    end
                end
                else if (window_last) begin
                    // Close, this event included. block_synced is pulsed WITH
                    // block_pulse. It must not be an internal state that the
                    // reporter samples one clock later, or every packet reports
                    // zero and the host rejects the whole capture.
                    block_pulse  <= 1'b1;
                    block_synced <= 1'b1;
                    block_id     <= next_block_id;

                    attempts   <= c_att   + i_att;
                    clicks_all <= c_click + i_click;
                    lost_all   <= c_lost  + i_lost;
                    trials_z   <= c_tz + i_tz;
                    clicks_z   <= c_cz + i_cz;
                    errors_z   <= c_ez + i_ez;
                    trials_x   <= c_tx + i_tx;
                    clicks_x   <= c_cx + i_cx;
                    errors_x   <= c_ex + i_ex;

                    decoy_counts  <= class_add(c_dec, i_dec);
                    vacuum_counts <= class_add(c_vac, i_vac);

                    h_mean <= hsum_n[27:ATTEMPT_LOG2];
                    h_min  <= hmin_n;
                    h_max  <= hmax_n;

                    blk_cfg        <= c_cfg;
                    config_changed <= c_cfg_chg | cfg_mismatch;
                    next_block_id  <= next_block_id + 32'd1;

                    // Re-arm. The next attempt is window_first and re-opens
                    // the block through the !armed branch, which latches the
                    // configuration afresh. Staying armed here would compare
                    // the new block's first attempt against the OLD c_cfg, so a
                    // controller change landing between two blocks would flag
                    // the next, perfectly homogeneous, block as changed.
                    armed   <= 1'b0;
                    c_att   <= 18'd0; c_click <= 18'd0; c_lost <= 18'd0;
                    c_tz    <= 18'd0; c_cz    <= 18'd0; c_ez   <= 18'd0;
                    c_tx    <= 18'd0; c_cx    <= 18'd0; c_ex   <= 18'd0;
                    c_dec   <= 144'd0; c_vac  <= 144'd0;
                    c_hsum  <= 28'd0; c_hmin  <= 12'hFFF; c_hmax <= 12'd0;
                    c_cfg_chg <= 1'b0;
                end
                else begin
                    c_att   <= c_att   + i_att;
                    c_click <= c_click + i_click;
                    c_lost  <= c_lost  + i_lost;
                    c_tz    <= c_tz + i_tz;
                    c_cz    <= c_cz + i_cz;
                    c_ez    <= c_ez + i_ez;
                    c_tx    <= c_tx + i_tx;
                    c_cx    <= c_cx + i_cx;
                    c_ex    <= c_ex + i_ex;
                    c_dec   <= class_add(c_dec, i_dec);
                    c_vac   <= class_add(c_vac, i_vac);
                    c_hsum  <= hsum_n;
                    c_hmin  <= hmin_n;
                    c_hmax  <= hmax_n;
                    if (cfg_mismatch) c_cfg_chg <= 1'b1;
                end
            end
        end
    end

endmodule
