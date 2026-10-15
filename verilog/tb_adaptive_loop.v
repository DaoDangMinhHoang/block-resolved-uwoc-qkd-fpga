// ============================================================
// TESTBENCH — closed loop: uwoc_channel + channel_monitor + adaptive_controller
// ============================================================
// Validates STEP 3. Skips the OOK TX/RX and the top_module FSM, wiring the
// qubit-level events from the channel straight into the monitor — enough to
// exercise the control loop.
//
// TEST A — SNR normalisation: at the nominal operating point (mu_level = 8 =
//          mu_ref) snr_level must be ~128 REGARDLESS of water type / distance.
//          That is the whole point of the normalisation: remove the static
//          exp(-c*d) attenuation.
// TEST B — window_valid gate: long distance -> too few clicks -> window_valid = 0
//          -> the controller MUST HOLD its mode instead of chasing noise.
// TEST C — operating range: the mode must degrade as distance grows.
// TEST D — wavelength hill climbing: harbor must converge to red (lam=2),
//          clear ocean must stay on blue (lam=0).
// TEST E — [v14] LAMBDA LOCK-UP GUARD. Within ONE measurement point (no reset)
//          the hill climber still probes candidate wavelengths periodically. In
//          clear ocean at 15 m the 650 nm candidate drops P_click by 27x -> the
//          window stops being valid -> the lambda state machine (which sits
//          inside the window_valid branch) STALLS FOREVER and keeps the bad
//          wavelength. This is exactly the bug that wrecked the overnight
//          measurement sessions: the last 4 points of the 2026-08-09 session
//          measured P_click = 2.03/1.42/1.37/1.35e-5, matching the 650 nm
//          prediction (1.94/1.41/1.32/1.30e-5) precisely.
//          TEST E requires most windows to stay VALID across 60 windows.
// ============================================================
`timescale 1ns / 1ps

module tb_adaptive_loop;

    // Shortened window (2^15) to keep the simulation fast; NEXP_LOG2 stays 16
    // (the window size used when the LUT was generated) and the monitor
    // compensates for the difference with a bit shift — so the ROM does NOT
    // have to be regenerated just for simulation.
    localparam AW   = 15;
    localparam NEXP = 16;
    localparam integer WIN = (1 << AW);

    reg clk = 1'b0, rst_n = 1'b0;
    always #10 clk = ~clk;

    // ---- configuration ----
    reg  [1:0] water     = 2'd0;
    reg  [3:0] dist      = 4'd0;
    reg  [2:0] turb      = 3'd1;
    reg        adaptive  = 1'b0;
    reg  [3:0] man_power = 4'd8;

    // [v14] Host 0x01 pulse: fpga_collect.py calls Link.configure() before EVERY
    // measurement point, and the trailing 0x01 byte raises pc_reset_req -> which
    // drives both channel_monitor.clear and adaptive_controller.cfg_rst. The TB
    // must model this exactly, otherwise it is testing a scenario that never
    // happens.
    reg        cfg_rst   = 1'b0;

    // ---- channel ----
    reg         sample_en = 1'b0;
    wire        click, no_click, err_inject;
    wire [15:0] nexp_inv;
    wire [2:0]  cur_level;
    wire [11:0] h_s, h_o, h_f;

    // ---- controller -> channel ----
    wire [3:0] adapt_power;
    wire [7:0] adapt_basis;
    wire [1:0] adapt_lambda;      // lambda being transmitted (may be a candidate under probe)
    wire [1:0] adapt_lam_best;    // CONVERGED lambda
    wire [1:0] adapt_mode;
    wire       tx_allowed;
    wire [7:0] key_rate, stale;
    wire [23:0] adapt_slot;
    wire [7:0]  adapt_gap;

    wire [3:0] eff_power  = adaptive ? adapt_power  : man_power;
    wire [1:0] eff_lambda = adaptive ? adapt_lambda : 2'd0;

    uwoc_channel #(.COH_LOG2(12)) chan (
        .clk(clk), .rst_n(rst_n),
        .signal_in(1'b0), .signal_out(),
        .slot_width(24'd20),
        .chan_enable(1'b1),
        .water_type(water), .dist_idx(dist), .turb_level(turb),
        .lambda_sel(eff_lambda), .mu_level(eff_power),
        .p_noise(24'd218), .dynamic_enable(1'b0),
        .mu_gain_log2(2'd0), .decoy_enable(1'b0),
        .att_class(2'd0), .decoy_shift(3'd0),
        .current_level(cur_level),
        .click(click), .no_click(no_click), .err_inject(err_inject),
        .click_count(), .lost_count(), .flip_count(),
        .h_s_out(h_s), .h_o_out(h_o), .h_f_out(h_f),
        .nexp_inv(nexp_inv),
        .sample_en(sample_en),
        // coh_sel = 0 → legacy clock timer, so the controller regression tests keep
        // the fading cadence they were tuned against. The qubit-counted mode has a
        // dedicated check in tb_uwoc_channel.v (check 9).
        .coh_sel(4'd0),
        .stats_rst(1'b0)
    );

    // ---- basis-selection model (done inside the TB) ----
    reg [15:0] lfsr_b = 16'hACE1;
    always @(posedge clk)
        lfsr_b <= {lfsr_b[14:0], lfsr_b[15]^lfsr_b[13]^lfsr_b[12]^lfsr_b[10]};

    wire [7:0] basis_prob = adaptive ? adapt_basis : 8'd128;
    wire a_z = (lfsr_b[7:0]   < basis_prob);
    wire b_z = (lfsr_b[15:8]  < basis_prob);
    wire basis_match = (a_z == b_z);

    // ---- monitor ----
    wire [7:0]  qber, snr, photon, sifted, loss, jitter;
    wire [15:0] photon_cnt;
    wire        win_valid, win_pulse;

    channel_monitor #(.ATTEMPT_LOG2(AW), .NEXP_LOG2(NEXP), .MIN_SIFT(16)) mon (
        .clk(clk), .rst_n(rst_n),
        .evt_qubit_done(click),
        .evt_qubit_lost(no_click),
        .evt_basis_match(basis_match),
        .evt_data_error(err_inject & basis_match),
        .nexp_inv(nexp_inv),
        .enable(1'b1), .clear(cfg_rst),
        .qber(qber), .snr_level(snr),
        .photon_rate(photon), .photon_count(photon_cnt),
        .sifted_rate(sifted),
        .loss_rate(loss), .qber_jitter(jitter),
        .window_valid(win_valid), .window_pulse(win_pulse)
    );

    adaptive_controller ctrl (
        .clk(clk), .rst_n(rst_n),
        .qber(qber), .snr_level(snr), .photon_rate(photon),
        .photon_count(photon_cnt),
        .qber_jitter(jitter), .loss_rate(loss),
        .window_valid(win_valid), .window_pulse(win_pulse),
        .adaptive_enable(adaptive),
        .cfg_rst(cfg_rst),
        .manual_power(man_power), .manual_basis_prob(8'd128),
        .manual_slot_width(24'd250000), .manual_lambda(2'd0),
        .power_level(adapt_power), .basis_prob_z(adapt_basis),
        .slot_width_out(adapt_slot), .rep_gap(adapt_gap),
        .lambda_sel(adapt_lambda), .lambda_best(adapt_lam_best),
        .mode(adapt_mode),
        .tx_allowed(tx_allowed), .key_rate_est(key_rate),
        .stale_windows(stale)
    );

    integer i, fails;

    // window_pulse is high for EXACTLY 1 cycle and falls in the gap between two
    // pulse_qubit samples -> it must be counted in a dedicated always block; it
    // cannot be polled from inside a task.
    integer win_count, valid_count;
    initial begin win_count = 0; valid_count = 0; end
    always @(posedge clk) if (win_pulse) begin
        win_count = win_count + 1;
        if (win_valid) valid_count = valid_count + 1;
    end

    // Free-running qubit generator, independent of the TB control flow.
    initial begin
        sample_en = 1'b0;
        @(posedge rst_n);
        forever begin
            @(negedge clk); sample_en = 1'b1;
            @(negedge clk); sample_en = 1'b0;
        end
    end

    // Wait until n_win windows have elapsed
    task run_windows;
        input integer n_win;
        integer target;
    begin
        target = win_count + n_win;
        wait (win_count >= target);
        @(negedge clk);
    end
    endtask

    // [v14] Start a new MEASUREMENT POINT exactly the way the host's
    // Link.configure() does: load the configuration, then emit the 0x01 pulse.
    // Without that pulse the controller state (mode, lambda, stale) carries over
    // from the previous point — and one bad point then poisons the whole session.
    task new_point;
        input [1:0] w;
        input [3:0] d;
        input [2:0] lv;
    begin
        water = w; dist = d; turb = lv;
        @(negedge clk); cfg_rst = 1'b1;
        @(negedge clk); cfg_rst = 1'b0;
    end
    endtask

    task show;
        input [8*24:1] tag;
    begin
        $display("  %0s w=%0d d=%2d lv=%0d | snr=%3d qber=%3d(%.1f%%) jit=%3d clicks=%5d sift=%3d valid=%b | mode=%0d lam=%0d(best=%0d) pow=%2d stale=%0d",
                 tag, water, dist, turb, snr, qber, qber/2.0, jitter,
                 photon_cnt, sifted, win_valid, adapt_mode, adapt_lambda,
                 adapt_lam_best, adapt_power, stale);
    end
    endtask

    initial begin
        fails = 0;
        $display("");
        $display("================================================================");
        $display("  TB vong kin: uwoc_channel + channel_monitor + adaptive_controller");
        $display("  Cua so = %0d attempt (ATTEMPT_LOG2=%0d)", WIN, AW);
        $display("================================================================");

        repeat (5) @(negedge clk); rst_n = 1'b1; repeat (20) @(negedge clk);

        // ================= TEST A =================
        $display("");
        $display("[TEST A] SNR chuan hoa ~128 tai diem danh dinh (mu_level=8)");
        $display("  Muc dich: loai bo suy hao tinh exp(-c*d) khoi chi so SNR.");
        adaptive = 1'b0; man_power = 4'd8;
        begin : testA
            integer nA, okA;
            okA = 0; nA = 0;
            // clear ocean 5 m / 15 m / 25 m ; coastal 2 m / 6 m ; harbor 0.5 m
            for (i = 0; i < 6; i = i + 1) begin
                case (i)
                  0: begin water=2'd0; dist=4'd0; end
                  1: begin water=2'd0; dist=4'd2; end
                  2: begin water=2'd0; dist=4'd4; end
                  3: begin water=2'd1; dist=4'd0; end
                  4: begin water=2'd1; dist=4'd2; end
                  5: begin water=2'd2; dist=4'd0; end
                endcase
                turb = 3'd1;
                run_windows(3);
                show("A:");
                nA = nA + 1;
                // accept 96..176 (128 +/- 25%): finite sampling + ROM rounding
                if (snr >= 96 && snr <= 176) okA = okA + 1;
                else $display("     *** snr=%0d nam ngoai [96,176]", snr);
            end
            if (okA >= 5)
                $display("  [OK] %0d/%0d diem co snr trong [96,176] quanh 128", okA, nA);
            else begin
                $display("  *** FAIL: chi %0d/%0d diem dat", okA, nA);
                fails = fails + 1;
            end
        end

        // ================= TEST B =================
        $display("");
        $display("[TEST B] Cong window_valid: cu ly xa -> khong du click");
        adaptive = 1'b1;
        water = 2'd0; dist = 4'd10; turb = 3'd1;    // clear ocean 55 m: dead link
        run_windows(6);
        show("B:");
        if (!win_valid && stale > 0) begin
            $display("  [OK] window_valid=0 va stale=%0d -> controller khong hanh dong theo nhieu",
                     stale);
        end else begin
            $display("  *** FAIL: ky vong window_valid=0 o cu ly da chet link");
            fails = fails + 1;
        end
        // after many silent windows -> must fall back to PAUSE
        run_windows(10);
        show("B2:");
        if (adapt_mode == 2'd3)
            $display("  [OK] mat link keo dai -> mode=PAUSE(3), tx_allowed=%b", tx_allowed);
        else begin
            $display("  *** FAIL: mat link nhung mode=%0d (ky vong 3)", adapt_mode);
            fails = fails + 1;
        end

        // ================= TEST C =================
        $display("");
        $display("[TEST C] Pham vi hoat dong: mode xau dan khi cu ly tang");
        adaptive = 1'b1; turb = 3'd1; water = 2'd0;
        begin : testC
            integer m0, m2, m4, m6;
            dist = 4'd0; run_windows(12); show("C:"); m0 = adapt_mode;
            dist = 4'd2; run_windows(12); show("C:"); m2 = adapt_mode;
            dist = 4'd4; run_windows(12); show("C:"); m4 = adapt_mode;
            dist = 4'd6; run_windows(12); show("C:"); m6 = adapt_mode;
            // STRICT requirement: at long range (35 m, beyond the model's
            // d_max ~28 m) the controller must move DEFINITIVELY to CONSERVATIVE
            // or worse; it must not sit at MODERATE while it never collects
            // enough samples to conclude anything.
            if (m6 > m0 && m6 >= 2)
                $display("  [OK] mode(5m)=%0d -> mode(35m)=%0d: xau dan va dat >= CONSERVATIVE",
                         m0, m6);
            else begin
                $display("  *** FAIL: mode khong xau dan dung (%0d %0d %0d %0d)",
                         m0, m2, m4, m6);
                fails = fails + 1;
            end
        end

        // ================= TEST D =================
        $display("");
        $display("[TEST D] Leo doi buoc song (0=450 blue, 1=532 green, 2=650 red)");
        $display("  Ky vong Kebapci'23 §I: clear ocean -> blue; nuoc duc -> red.");
        adaptive = 1'b1;

        new_point(2'd2, 4'd0, 3'd1);                // harbor 0.5 m
        run_windows(60);
        show("D harbor:");
        // Check lambda_best (the CONVERGED result), not lambda_sel — while a
        // probe is running lambda_sel temporarily holds the candidate lambda.
        if (adapt_lam_best == 2'd2)
            $display("  [OK] harbor hoi tu ve lam_best=2 (650 nm red)");
        else begin
            $display("  *** FAIL: harbor hoi tu ve lam_best=%0d (ky vong 2)",
                     adapt_lam_best);
            fails = fails + 1;
        end

        // NEXT measurement point -> the host emits 0x01. Skip this step and
        // lam_best = 2 (correct for harbor) carries over into clear ocean, where
        // 650 nm kills the link immediately and no valid window is left for the
        // hill climber to recover from.
        new_point(2'd0, 4'd2, 3'd1);                // clear ocean 15 m
        run_windows(60);
        show("D clear:");
        if (adapt_lam_best == 2'd0)
            $display("  [OK] clear ocean hoi tu ve lam_best=0 (450 nm blue)");
        else begin
            $display("  *** FAIL: clear ocean hoi tu ve lam_best=%0d (ky vong 0)",
                     adapt_lam_best);
            fails = fails + 1;
        end

        // ================= TEST E =================
        $display("");
        $display("[TEST E] [v14] Do ung vien lam chet link -> phai HUY, khong duoc ket");
        $display("  clear ocean 15 m: 450 nm cho ~57 bit sang/cua so (hop le),");
        $display("  650 nm chi con ~2 (khong hop le). Bo leo doi van phai chay tiep.");
        adaptive = 1'b1;
        new_point(2'd0, 4'd2, 3'd1);
        begin : testE
            integer v0, w0, nval, nwin;
            v0 = valid_count; w0 = win_count;
            run_windows(60);
            nval = valid_count - v0;
            nwin = win_count - w0;
            show("E:");
            $display("  cua so hop le: %0d/%0d", nval, nwin);
            // With the old bug: the first window that lands on 650 nm becomes
            // invalid FOREVER -> nval stops at the first few tens of percent.
            // With the fix: only exactly 1 window is lost per 650 nm probe cycle.
            if (nval * 4 >= nwin * 3 && adapt_lam_best == 2'd0)
                $display("  [OK] %0d/%0d cua so hop le va lam_best=0 -> khong bi ket",
                         nval, nwin);
            else begin
                $display("  *** FAIL: chi %0d/%0d cua so hop le, lam_best=%0d (ky vong >=75%% va 0)",
                         nval, nwin, adapt_lam_best);
                fails = fails + 1;
            end
        end

        $display("");
        $display("================================================================");
        if (fails == 0) $display("  KET QUA: TAT CA TEST DAT");
        else            $display("  KET QUA: %0d TEST HONG", fails);
        $display("================================================================");
        $display("");
        $finish;
    end

endmodule
