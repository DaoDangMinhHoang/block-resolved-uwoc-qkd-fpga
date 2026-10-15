// ============================================================
// ADAPTIVE CONTROLLER v2 — QKD parameter optimisation for the UWOC channel
// ============================================================
//
// THIS IS THE CORE CONTRIBUTION OF THE PAPER.
//
// CHANGES FROM v1 (the atmospheric-FSO version) — every item follows from
// a quantitative conclusion of python/uwoc_channel_model.py:
//
//  [1] window_valid GATE.  Underwater a window may collect only a few clicks;
//      the QBER is then pure shot noise. v1 would switch mode on that noise.
//      Now, if !window_valid the mode is HELD unchanged (CHECK 7).
//
//  [2] qber_jitter AS THE TURBULENCE INDICATOR.  CHECK 8: across 5 turbulence
//      levels the MEAN QBER barely moves (70% variation) while the SPREAD
//      between windows changes 252%. Looking at mean QBER alone, the
//      controller would be BLIND to underwater turbulence. The real risk is
//      OUTAGE (a stray window past 11%), so high jitter ⇒ pre-emptive downgrade.
//
//  [3] SNR THRESHOLDS ON THE NORMALISED SCALE.  channel_monitor v3 outputs
//      snr_level = 128·(actual / expected clicks), i.e. the FADING MARGIN with
//      static loss exp(−c·d) removed. The old thresholds (180/100/40 on the
//      "qubits received" scale) no longer mean anything.
//
//  [4] NEW KNOB: WAVELENGTH λ.  This is the adaptation axis SPECIFIC to water
//      and the only knob acting directly on c(λ) — the dominant attenuation
//      term. Implemented by hill-climbing on the measured click rate.
//      Expectation (Kebapci'23 §I): clear ocean → blue-green; turbid water →
//      shift toward red-yellow. RTL check: TEST 1b of tb_uwoc_channel.v.
//
//  [5] POWER CEILING.  v1 set CON_POWER = 15 (max μ). Raising μ does lower
//      QBER (QBER → e_pol once n̄ ≫ Y₀) but it SIMULTANEOUSLY raises the
//      multi-photon pulse fraction ⇒ opening the door to PNS attacks. The huge
//      underwater loss hides PNS even better, so μ is capped at MU_CAP.
//
// 4 MODES (semantics unchanged from v1)
//   AGGRESSIVE (00) good channel → basis 50/50, narrow slot, low μ, max rate
//   MODERATE   (01) average      → basis 60/40, balanced
//   CONSERVATIVE(10) bad channel → basis 80/20, wide slot, high μ, reliability
//   PAUSE      (11) link lost    → stop transmitting, protect the key
// ============================================================

module adaptive_controller #(
    // Transmit intensity ceiling (anti-PNS). 12/15 ≈ 1.5× nominal μ.
    parameter [3:0] MU_CAP = 4'd12,
    // Settling windows before each reference measurement of the wavelength prober.
    // One probe cycle = LAM_HOLD + 2·LAM_ACC + 1 windows; converging from the
    // initial λ to the optimum takes at most 3 cycles (there are 3 λ choices).
    parameter [3:0] LAM_HOLD = 4'd2
) (
    input  wire        clk,
    input  wire        rst_n,

    // ---- Measurements from channel_monitor v3 ----
    input  wire [7:0]  qber,           // 0–200, 0.5%/unit
    input  wire [7:0]  snr_level,      // NORMALISED: 128 = nominal
    input  wire [7:0]  photon_rate,    // clicks/window (saturates at 255)
    input  wire [15:0] photon_count,   // clicks/window, NOT saturating — used by
                                       // the λ hill-climber (see its section)
    input  wire [7:0]  qber_jitter,    // EWMA |ΔQBER| — turbulence indicator
    input  wire [7:0]  loss_rate,      // 255 = everything lost
    input  wire        window_valid,   // estimate trustworthy?
    input  wire        window_pulse,

    input  wire        adaptive_enable,

    // [v14] Pulse from the host's 0x01 "reset statistics" command. Re-arms the
    // WHOLE controller — mode, λ, stale counter — so one measurement point cannot
    // hand its state to the next.
    //
    // ⚠ WITHOUT THIS THE OVERNIGHT SESSIONS WERE UNRECOVERABLE. The host runs 10
    //   points back to back, sending 0x01 between them; that pulse cleared
    //   channel_monitor and uwoc_channel but NOT this module, whose only reset was
    //   KEY[3] on the board. So a controller that got stuck at the wrong λ during
    //   point 3 stayed stuck for points 4…10 as well, and the whole session past
    //   the first failure was collected at the noise floor.
    input  wire        cfg_rst,

    // Manual override (when adaptive_enable = 0)
    input  wire [3:0]  manual_power,
    input  wire [7:0]  manual_basis_prob,
    input  wire [23:0] manual_slot_width,
    input  wire [1:0]  manual_lambda,

    // ---- Control outputs ----
    output reg  [3:0]  power_level,    // → uwoc_channel's mu_level
    output reg  [7:0]  basis_prob_z,   // 128=50%, 154=60%, 204=80%
    output reg  [23:0] slot_width_out,
    output reg  [7:0]  rep_gap,
    output reg  [1:0]  lambda_sel,     // λ IN USE — while probing this is the
                                       // candidate, NOT the converged result
    output wire [1:0]  lambda_best,    // converged λ — for reporting/evaluation
    output reg  [1:0]  mode,
    output reg         tx_allowed,
    output reg  [7:0]  key_rate_est,
    output reg  [7:0]  stale_windows   // consecutive windows with too few samples
);

    // ============================
    // Thresholds
    // ============================
    // QBER, 0–200 scale (0.5%/unit)
    localparam [7:0] QBER_GOOD = 8'd8;    //  4%
    localparam [7:0] QBER_WARN = 8'd16;   //  8%
    localparam [7:0] QBER_BAD  = 8'd22;   // 11% — BB84 security limit
    localparam [7:0] QBER_CRIT = 8'd30;   // 15%

    // SNR on the NORMALISED scale (128 = nominal click rate)
    localparam [7:0] SNR_HIGH = 8'd160;   // 1.25× nominal
    localparam [7:0] SNR_MED  = 8'd96;    // 0.75×
    localparam [7:0] SNR_LOW  = 8'd40;    // 0.31× — channel almost dead

    // Jitter (same scale as QBER, 0.5%/unit).
    //
    // ⚠ THE THRESHOLD MUST SIT ABOVE THE SHOT-NOISE FLOOR. With ~90 sifted
    //   samples and QBER ~2%, the binomial spread is √(0.02·0.98/90) ≈ 1.5%,
    //   so |ΔQBER| between two windows already averages ≈ 3–4 units EVEN ON A
    //   COMPLETELY STATIC CHANNEL. The first trial (JIT_HIGH = 6) sat below
    //   that floor → spurious CONSERVATIVE drops; caught by tb_adaptive_loop.v.
    //
    //   Calibrated against CHECK 8 (uwoc_channel_model.py), same sample size:
    //     weak turbulence   → std(QBER) ≈  1.74%  →  ~3.5 units
    //     SEVERE turbulence → std(QBER) ≈ 11.32%  →  ~23  units
    //   ⇒ put the thresholds in between, well clear of the shot-noise floor.
    localparam [7:0] JIT_LOW  = 8'd6;     // genuinely static channel
    localparam [7:0] JIT_HIGH = 8'd16;    // real turbulence → pre-emptive downgrade

    // Consecutive windows with NO clicks at all before declaring the link lost.
    localparam [7:0] DEAD_WINDOWS = 8'd8;

    // [v14] "No clicks at all" CANNOT be tested as photon_rate == 0. Y₀ (dark +
    // background) is 1.30e-5 per attempt, so a 65 536-attempt window collects
    // 0.85 noise clicks on average even with the transmitter silent — the test is
    // false about 57% of the time on a link that is completely dead. Compare
    // against the noise floor instead: P(Poisson(0.85) > 3) < 1%, while the
    // weakest range the sweep visits still delivers ~6 clicks/window.
    localparam [7:0] DEAD_CLICKS = 8'd3;

    // ============================
    // Per-mode parameters
    // ============================
    localparam [3:0]  AGG_POWER = 4'd6;
    localparam [7:0]  AGG_BASIS = 8'd128;      // 50% Z
    localparam [23:0] AGG_SLOT  = 24'd250000;  // 5 ms
    localparam [7:0]  AGG_GAP   = 8'd1;

    localparam [3:0]  MOD_POWER = 4'd9;
    localparam [7:0]  MOD_BASIS = 8'd154;      // 60% Z
    localparam [23:0] MOD_SLOT  = 24'd500000;  // 10 ms
    localparam [7:0]  MOD_GAP   = 8'd2;

    localparam [3:0]  CON_POWER = MU_CAP;      // capped, anti-PNS
    localparam [7:0]  CON_BASIS = 8'd204;      // 80% Z
    localparam [23:0] CON_SLOT  = 24'd2500000; // 50 ms
    localparam [7:0]  CON_GAP   = 8'd4;

    localparam [7:0]  PAU_GAP   = 8'd16;

    // ============================
    // Mode decision
    // ============================
    reg [2:0] upgrade_counter;

    // Purely combinational: desired mode from the current measurements.
    reg [1:0] want_mode;
    always @(*) begin
        if (qber >= QBER_CRIT || snr_level <= SNR_LOW)
            want_mode = 2'b11;                                   // PAUSE
        else if (qber >= QBER_WARN || snr_level < SNR_MED
                 || qber_jitter >= JIT_HIGH)
            want_mode = 2'b10;                                   // CONSERVATIVE
        else if (qber < QBER_GOOD && snr_level >= SNR_HIGH
                 && qber_jitter < JIT_LOW)
            want_mode = 2'b00;                                   // AGGRESSIVE
        else
            want_mode = 2'b01;                                   // MODERATE
    end

    // ============================
    // Wavelength probing by hill-climbing
    // ============================
    // Cycle: SETTLE (settle on the best λ) → BASE (accumulate the reference)
    //        → WARM (switch λ, drop 1 window) → MEAS (accumulate candidate) → decide
    //
    // ⚠ TWO MANDATORY ADJUSTMENTS, found by tb_adaptive_loop.v TEST D:
    //
    //  (a) USE the 16-bit photon_count, NOT the 8-bit photon_rate.
    //      photon_rate saturates at 255; at short range every λ hits the cap, so
    //      comparing between wavelengths becomes meaningless.
    //
    //  (b) ACCUMULATE OVER SEVERAL WINDOWS. The true gap between two adjacent λ
    //      can be only ~10% (e.g. 532 nm vs 650 nm in harbour water: 9.9%), while
    //      the standard deviation of a ~220-click window is 1/√220 ≈ 6.7%. With a
    //      1/8 = 12.5% acceptance margin on ONE window, the hill-climber CANNOT
    //      resolve the real gap and stalls on the wrong λ. Averaging over
    //      LAM_ACC = 4 windows brings the spread down to ~3.35%, letting the
    //      margin drop to 1/16 = 6.25% (≈3σ) without chasing measurement noise.
    localparam [2:0] LAM_SETTLE = 3'd0,
                     LAM_BASE   = 3'd1,
                     LAM_WARM   = 3'd2,
                     LAM_MEAS   = 3'd3;

    localparam [3:0] LAM_ACC = 4'd4;      // windows accumulated per measurement

    reg [2:0]  lam_state;
    reg [1:0]  lam_best;         // λ currently believed best
    reg [1:0]  lam_home;         // [v14] last λ that actually produced a VALID
                                 // window while it was the converged choice —
                                 // the value to fall back to when a probe kills
                                 // the link. Updated only in SETTLE/BASE, where
                                 // lambda_sel == lam_best, so a probe candidate
                                 // can never be mistaken for the converged λ.
    reg [1:0]  lam_cand;         // λ under test
    reg [19:0] lam_base_acc;     // total clicks at lam_best
    reg [19:0] lam_cand_acc;     // total clicks at lam_cand
    reg [3:0]  lam_timer;
    reg        lam_dir;          // 0 = try λ+1, 1 = try λ−1

    wire [1:0] lam_up   = (lam_best == 2'd2) ? 2'd0 : (lam_best + 2'd1);
    wire [1:0] lam_down = (lam_best == 2'd0) ? 2'd2 : (lam_best - 2'd1);

    assign lambda_best = lam_best;

    // ============================
    // Main sequential block
    // ============================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            mode            <= 2'b01;
            power_level     <= MOD_POWER;
            basis_prob_z    <= MOD_BASIS;
            slot_width_out  <= MOD_SLOT;
            rep_gap         <= MOD_GAP;
            lambda_sel      <= 2'd0;
            tx_allowed      <= 1'b1;
            key_rate_est    <= 8'd0;
            upgrade_counter <= 3'd0;
            stale_windows   <= 8'd0;
            lam_state       <= LAM_SETTLE;
            lam_best        <= 2'd0;
            lam_home        <= 2'd0;
            lam_cand        <= 2'd0;
            lam_base_acc    <= 20'd0;
            lam_cand_acc    <= 20'd0;
            lam_timer       <= 4'd0;
            lam_dir         <= 1'b0;
        end else if (cfg_rst) begin
            // ---- [v14] Host asked for a fresh measurement point ----
            // Same values as the asynchronous reset. It MUST be its own branch and
            // may NOT be folded into `if (!rst_n || cfg_rst)`: Quartus only accepts
            // signals present in the sensitivity list inside an async-reset branch
            // (Error 10200), exactly as documented for `clear` in channel_monitor.v.
            mode            <= 2'b01;
            power_level     <= MOD_POWER;
            basis_prob_z    <= MOD_BASIS;
            slot_width_out  <= MOD_SLOT;
            rep_gap         <= MOD_GAP;
            lambda_sel      <= 2'd0;
            tx_allowed      <= 1'b1;
            key_rate_est    <= 8'd0;
            upgrade_counter <= 3'd0;
            stale_windows   <= 8'd0;
            lam_state       <= LAM_SETTLE;
            lam_best        <= 2'd0;
            lam_home        <= 2'd0;
            lam_cand        <= 2'd0;
            lam_base_acc    <= 20'd0;
            lam_cand_acc    <= 20'd0;
            lam_timer       <= 4'd0;
            lam_dir         <= 1'b0;
        end else if (!adaptive_enable) begin
            // ---- Manual mode ----
            power_level    <= manual_power;
            basis_prob_z   <= manual_basis_prob;
            slot_width_out <= manual_slot_width;
            lambda_sel     <= manual_lambda;
            rep_gap        <= 8'd1;
            tx_allowed     <= 1'b1;
            mode           <= 2'b01;
            stale_windows  <= 8'd0;
        end else if (window_pulse) begin

            // ------------------------------------------------
            // (A) Window with too few samples → do NOT act on noise
            // ------------------------------------------------
            if (!window_valid) begin
                if (stale_windows != 8'd255)
                    stale_windows <= stale_windows + 1'b1;

                // ---- [v14] ABORT AN IN-FLIGHT λ PROBE ----
                //
                // THIS IS THE FAILURE THAT DESTROYED THE OVERNIGHT ADAPTIVE RUNS.
                // The hill-climber APPLIES the candidate (LAM_WARM sets
                // lambda_sel <= lam_cand) and only reverts it at the end of
                // LAM_MEAS — but the whole λ FSM lives inside the window_valid
                // branch. Pick a candidate bad enough to kill the link and the
                // window stops being valid, so the FSM never reaches its own
                // revert and the BAD λ STAYS APPLIED FOREVER.
                //
                // In clear ocean 650 nm is 27× worse than 450 nm at 15 m and 92×
                // at 25 m, so the probe drops P_click to the Y₀ floor. Measured
                // against the ROM, the four dead points of the 2026-08-09 session
                // sat at 2.03/1.42/1.37/1.35e-5 where λ = 650 nm predicts
                // 1.94/1.41/1.32/1.30e-5 — the link was running on the probe
                // wavelength, not on 450 nm.
                //
                // v1's answer was the escape hatch below, but it only ran when
                // photon_rate was exactly 0, which the noise floor makes rare
                // (see DEAD_CLICKS). Aborting the probe directly is the real fix:
                // one invalid window and λ is back on the known-good value.
                if (lam_state != LAM_SETTLE) begin
                    lambda_sel <= lam_home;
                    lam_best   <= lam_home;
                    lam_state  <= LAM_SETTLE;
                    lam_timer  <= 4'd0;
                    lam_dir    <= ~lam_dir;   // try the OTHER direction next time
                end

                // ⚠ Do NOT use loss_rate as the "link lost" indicator underwater.
                //   loss_rate = fraction of attempts without a click, and underwater
                //   QKD has P_click ~ 10⁻³ ⇒ loss_rate SATURATES AT 255 even on a
                //   perfectly healthy link. The condition (loss_rate ≥ 250) is
                //   therefore ALWAYS TRUE and once pushed a healthy 15 m clear-ocean
                //   link into PAUSE — a bug found by tb_adaptive_loop.v. The correct
                //   indicator is a click count AT THE NOISE FLOOR.
                if (photon_rate <= DEAD_CLICKS) begin
                    // [v14] The old escape hatch ROTATED lam_best here, which at
                    // the far points (where the window is legitimately
                    // sample-starved on the correct λ) walked the converged
                    // wavelength off to 532 then 650 and left it there — a second,
                    // independent way to poison the rest of the session. λ now only
                    // ever falls BACK to the last value that produced a valid
                    // window; the probe abort above is what handles a bad
                    // candidate, so nothing needs to hunt blindly any more.
                    if (stale_windows >= DEAD_WINDOWS) begin
                        lambda_sel <= lam_home;
                        lam_best   <= lam_home;
                        mode       <= 2'b11;
                        tx_allowed <= 1'b0;
                        rep_gap    <= PAU_GAP;
                    end
                end else begin
                    // Clicks, but too few samples: the link is alive, just sparse.
                    tx_allowed <= 1'b1;

                    // THE ESTIMATOR IS SAMPLE-STARVED — the right response is to
                    // make samples arrive faster, not to sit still:
                    //   • a strongly biased basis (p_z = 80%) lifts sifting
                    //     efficiency 0.50 → 0.68 = +36% sifted qubits/window;
                    //   • a higher μ (CON_POWER) directly raises the click rate.
                    // Without this branch the controller sticks in MODERATE at
                    // long range forever and never gathers enough samples to
                    // conclude — exactly the case tb_adaptive_loop.v TEST C builds.
                    if (stale_windows >= 8'd4 && mode < 2'b10) begin
                        mode           <= 2'b10;          // CONSERVATIVE
                        power_level    <= CON_POWER;
                        basis_prob_z   <= CON_BASIS;
                        slot_width_out <= CON_SLOT;
                        rep_gap        <= CON_GAP;
                        upgrade_counter <= 3'd0;
                    end
                end
            end else begin
                stale_windows <= 8'd0;

                // [v14] Remember the λ that is carrying this valid window, but
                // ONLY while the applied λ is the converged one. In LAM_WARM /
                // LAM_MEAS lambda_sel is a candidate under test; recording it here
                // would let a candidate that merely survived one window be
                // promoted to "home" and defeat the whole point of the fallback.
                if (lam_state == LAM_SETTLE || lam_state == LAM_BASE)
                    lam_home <= lam_best;

                // Asymmetric hysteresis: downgrade IMMEDIATELY (safety first),
                // upgrade requires 3 consecutive good windows.
                if (want_mode > mode) begin
                    mode            <= want_mode;
                    upgrade_counter <= 3'd0;
                end else if (want_mode < mode) begin
                    upgrade_counter <= upgrade_counter + 1'b1;
                    if (upgrade_counter >= 3'd3) begin
                        mode            <= want_mode;
                        upgrade_counter <= 3'd0;
                    end
                end else begin
                    upgrade_counter <= 3'd0;
                end

                // -------- Apply per-mode parameters --------
                case (mode)
                    2'b00: begin
                        power_level    <= AGG_POWER;
                        basis_prob_z   <= AGG_BASIS;
                        slot_width_out <= AGG_SLOT;
                        rep_gap        <= AGG_GAP;
                        tx_allowed     <= 1'b1;
                    end
                    2'b01: begin
                        power_level    <= MOD_POWER;
                        basis_prob_z   <= MOD_BASIS;
                        slot_width_out <= MOD_SLOT;
                        rep_gap        <= MOD_GAP;
                        tx_allowed     <= 1'b1;
                    end
                    2'b10: begin
                        power_level    <= CON_POWER;
                        basis_prob_z   <= CON_BASIS;
                        slot_width_out <= CON_SLOT;
                        rep_gap        <= CON_GAP;
                        tx_allowed     <= 1'b1;
                    end
                    2'b11: begin
                        power_level    <= CON_POWER;
                        basis_prob_z   <= CON_BASIS;
                        slot_width_out <= CON_SLOT;
                        rep_gap        <= PAU_GAP;
                        tx_allowed     <= 1'b0;
                        // Keep a sparse probe going to detect channel recovery
                        if (snr_level > SNR_MED && qber < QBER_BAD)
                            tx_allowed <= 1'b1;
                    end
                endcase

                // -------- Key-rate estimate --------
                // key_rate ≈ sifted_rate · (1 − 2·H₂(QBER)); approximated by a
                // per-mode shift, good enough for display / logging.
                case (mode)
                    2'b00: key_rate_est <= photon_rate >> 1;
                    2'b01: key_rate_est <= photon_rate >> 2;
                    2'b10: key_rate_est <= photon_rate >> 3;
                    2'b11: key_rate_est <= 8'd0;
                endcase

                // ------------------------------------------------
                // (B) Hill-climbing on wavelength
                // ------------------------------------------------
                // Runs only when the measurement is trustworthy. The criterion
                // is CLICK RATE (photon_rate), which directly reflects c(λ)·d.
                case (lam_state)
                    LAM_SETTLE: begin
                        lambda_sel <= lam_best;
                        if (lam_timer >= LAM_HOLD) begin
                            lam_timer    <= 4'd0;
                            lam_base_acc <= 20'd0;
                            lam_state    <= LAM_BASE;
                        end else
                            lam_timer <= lam_timer + 1'b1;
                    end

                    LAM_BASE: begin
                        lam_base_acc <= lam_base_acc + {4'd0, photon_count};
                        if (lam_timer >= (LAM_ACC - 1'b1)) begin
                            lam_timer <= 4'd0;
                            lam_cand  <= lam_dir ? lam_down : lam_up;
                            lam_state <= LAM_WARM;
                        end else
                            lam_timer <= lam_timer + 1'b1;
                    end

                    LAM_WARM: begin
                        // Drop 1 window after switching λ so the measurement is
                        // not contaminated by data from the old wavelength.
                        lambda_sel   <= lam_cand;
                        lam_cand_acc <= 20'd0;
                        lam_state    <= LAM_MEAS;
                    end

                    LAM_MEAS: begin
                        lam_cand_acc <= lam_cand_acc + {4'd0, photon_count};
                        if (lam_timer >= (LAM_ACC - 1'b1)) begin
                            // Keep the candidate if it beats the 1/16 margin (~3σ
                            // over LAM_ACC windows) — sensitive enough for a real
                            // ~10% gap, tight enough not to chase measurement noise.
                            if ((lam_cand_acc + {4'd0, photon_count})
                                > (lam_base_acc + (lam_base_acc >> 4)))
                                lam_best <= lam_cand;
                            else
                                lambda_sel <= lam_best;   // revert
                            lam_dir   <= ~lam_dir;
                            lam_timer <= 4'd0;
                            lam_state <= LAM_SETTLE;
                        end else
                            lam_timer <= lam_timer + 1'b1;
                    end

                    default: lam_state <= LAM_SETTLE;
                endcase
            end
        end
    end

endmodule
