// ============================================================
// TESTBENCH — uwoc_channel.v
// ============================================================
// Cross-checks the RTL against the Python golden model (_rtl_exact_sim in uwoc_lut_gen.py).
//
// Run:
//   vlib work
//   vlog +incdir+. uwoc_channel.v tb_uwoc_channel.v
//   vsim -c -do "run -all; quit" tb_uwoc_channel
//
// TEST 1 — detection-stage statistics (§1–5 of uwoc_channel.v)
//          sweeps (water_type, dist_idx, turb_level), counting clicks / errors.
// TEST 2 — OOK signal path (§6): checks no_click → silence for the whole frame,
//          and err → flips ONLY the b0 data slot, leaving SYNC/basis untouched.
// TEST 3 — coherence time: h must stay FROZEN between two coh_ticks.
// ============================================================
`timescale 1ns / 1ps

module tb_uwoc_channel;

    // A small COH_LOG2 keeps simulation fast: fading is resampled every 256 cycles.
    // On real hardware the default of 18 is used (≈5.24 ms @50 MHz).
    localparam COH = 8;
    localparam integer N_QUBIT = 200000;

    reg clk = 1'b0, rst_n = 1'b0;
    always #10 clk = ~clk;              // 50 MHz

    reg         signal_in   = 1'b0;
    reg  [23:0] slot_width  = 24'd20;
    reg         chan_enable = 1'b1;
    reg  [1:0]  water_type  = 2'd0;
    reg  [3:0]  dist_idx    = 4'd0;
    reg  [2:0]  turb_level  = 3'd1;
    reg  [1:0]  lambda_sel  = 2'd0;     // 450 nm — khop mo hinh vang
    reg  [3:0]  mu_level    = 4'd8;
    reg  [23:0] p_noise     = 24'd218;  // Y0 = 1.30e-5 → round(1.3e-5·2^24)=218
    reg         dyn_en      = 1'b0;
    reg  [1:0]  mu_gain     = 2'd0;     // [rev B] signal gain 2^g
    reg         decoy_en    = 1'b0;     // [rev B] honour att_class
    reg  [1:0]  att_class   = 2'd0;     // [rev B] 0 signal, 1 decoy, 2 vacuum
    reg  [2:0]  decoy_shift = 3'd0;     // [rev B] decoy P_sig = P_sig >> shift
    reg         sample_en   = 1'b0;
    reg  [3:0]  coh_sel_r   = 4'd0;     // 0 = legacy clock timer (TEST 1–4)
    reg         stats_rst_r = 1'b0;

    wire        signal_out;
    wire [2:0]  cur_level;
    wire        click, no_click, err_inject;
    wire [15:0] click_cnt, lost_cnt, flip_cnt;
    wire [11:0] h_s, h_o, h_f;
    wire [15:0] nexp_inv;
    wire [23:0] qub_index;

    uwoc_channel #(.COH_LOG2(COH)) dut (
        .clk(clk), .rst_n(rst_n),
        .signal_in(signal_in), .signal_out(signal_out),
        .slot_width(slot_width),
        .chan_enable(chan_enable), .water_type(water_type),
        .dist_idx(dist_idx), .turb_level(turb_level),
        .lambda_sel(lambda_sel),
        .mu_level(mu_level), .p_noise(p_noise),
        .dynamic_enable(dyn_en),
        .mu_gain_log2(mu_gain), .decoy_enable(decoy_en),
        .att_class(att_class), .decoy_shift(decoy_shift),
        .current_level(cur_level),
        .click(click), .no_click(no_click), .err_inject(err_inject),
        .click_count(click_cnt), .lost_count(lost_cnt), .flip_count(flip_cnt),
        .h_s_out(h_s), .h_o_out(h_o), .h_f_out(h_f),
        .qub_index(qub_index),
        .nexp_inv(nexp_inv),
        .sample_en(sample_en),
        // Held at 0 (legacy clock timer) for TEST 1–4 so the golden-model comparison
        // measures exactly what it measured before [v12]; TEST 5 drives it.
        .coh_sel(coh_sel_r),
        .stats_rst(stats_rst_r)
    );

    integer n_click, n_err, n_lost, i;
    integer fails;

    // Emit one sample_en pulse (1 cycle) and capture the result.
    // TIMING: sample_pulse fires on the posedge BETWEEN the two negedges below;
    // the click/no_click/err_inject outputs are registered, so they are valid at
    // the second negedge and cleared on the next posedge — read them EXACTLY here.
    task pulse_qubit;
    begin
        @(negedge clk); sample_en = 1'b1;
        @(negedge clk); sample_en = 1'b0;
        #1;
        if (click)      n_click = n_click + 1;
        if (no_click)   n_lost  = n_lost  + 1;
        if (err_inject) n_err   = n_err   + 1;
    end
    endtask

    // exp_p, exp_q: gia tri vang tu mo hinh float (uwoc_channel_model.py,
    // KIEM CHUNG 4). Chi kiem tra chat khi so click du lon (>500), vi o cu ly
    // xa sai lech chu yeu la nhieu ban chu khong phai loi RTL.
    task run_case;
        input [1:0] wt;
        input [3:0] di;
        input [2:0] lv;
        input real  exp_p;      // P_click ky vong
        input real  exp_q;      // QBER ky vong (phan thap phan)
        real pclick, qber, dp, dq;
    begin
        water_type = wt; dist_idx = di; turb_level = lv;
        n_click = 0; n_err = 0; n_lost = 0;
        repeat (20) @(negedge clk);
        for (i = 0; i < N_QUBIT; i = i + 1) pulse_qubit;
        pclick = n_click * 1.0 / N_QUBIT;
        qber   = (n_click > 0) ? (n_err * 1.0 / n_click) : 0.5;
        dp = (pclick - exp_p) / exp_p * 100.0; if (dp < 0) dp = -dp;
        dq = (qber - exp_q) * 100.0;           if (dq < 0) dq = -dq;
        $display("  w=%0d d=%2d lv=%0d | eff=%0d | click=%6d lost=%6d err=%5d | P=%.3e (vang %.3e, %5.1f%%) QBER=%5.2f%% (vang %5.2f%%, %4.2f%%)",
                 wt, di, lv, cur_level, n_click, n_lost, n_err,
                 pclick, exp_p, dp, qber * 100.0, exp_q * 100.0, dq);
        if (n_click + n_lost != N_QUBIT) begin
            $display("    *** FAIL: click+lost=%0d != %0d", n_click + n_lost, N_QUBIT);
            fails = fails + 1;
        end
        if (n_click > 500) begin
            if (dp > 15.0) begin
                $display("    *** FAIL: P_click lech %.1f%% so voi mo hinh vang", dp);
                fails = fails + 1;
            end
            if (dq > 2.0) begin
                $display("    *** FAIL: QBER lech %.2f%% tuyet doi so voi vang", dq);
                fails = fails + 1;
            end
        end
    end
    endtask

    // ---- TEST 2: send one OOK frame and check the output waveform ----
    integer sync_hi, b1_seen, b0_seen, k;
    reg exp_b1, exp_b0;

    task send_frame;
        input b1, b0;
    begin
        exp_b1 = b1; exp_b0 = b0;
        sync_hi = 0; b1_seen = 0; b0_seen = 0;
        // SYNC = 1
        signal_in = 1'b1;
        for (k = 0; k < slot_width; k = k + 1) begin
            @(posedge clk); #1; if (signal_out) sync_hi = sync_hi + 1;
        end
        // slot b1 (basis)
        signal_in = b1;
        for (k = 0; k < slot_width; k = k + 1) begin
            @(posedge clk); #1; if (signal_out) b1_seen = b1_seen + 1;
        end
        // slot b0 (data)
        signal_in = b0;
        for (k = 0; k < slot_width; k = k + 1) begin
            @(posedge clk); #1; if (signal_out) b0_seen = b0_seen + 1;
        end
        // IDLE
        signal_in = 1'b0;
        repeat (slot_width) @(posedge clk);
    end
    endtask

    integer prev_hs, prev_ho, changes;

    initial begin
        fails = 0;
        $display("");
        $display("================================================================");
        $display("  TESTBENCH uwoc_channel.v   (COH_LOG2=%0d, %0d qubit/case)", COH, N_QUBIT);
        $display("================================================================");

        repeat (5) @(negedge clk);
        rst_n = 1'b1;
        repeat (10) @(negedge clk);

        // ---------- TEST 1 ----------
        $display("");
        $display("[TEST 1] Thong ke tang phat hien photon");
        $display("  water: 0=clear ocean, 1=coastal, 2=harbor");
        $display("  Doi chieu voi mo hinh vang Python (uwoc_lut_gen.py --verify)");
        $display("");
        //        water dist  lv   P_click vang   QBER vang
        run_case(2'd0, 4'd0, 3'd1, 1.0391e-2, 0.0164);  // clear ocean  5 m
        run_case(2'd0, 4'd0, 3'd4, 1.0391e-2, 0.0164);
        run_case(2'd0, 4'd2, 3'd1, 3.5108e-3, 0.0268);  // clear ocean 15 m
        run_case(2'd0, 4'd2, 3'd4, 3.5108e-3, 0.0268);
        run_case(2'd0, 4'd4, 3'd1, 1.1893e-3, 0.0369);  // clear ocean 25 m
        run_case(2'd0, 4'd6, 3'd1, 2.7805e-4, 0.0584);  // clear ocean 35 m
        run_case(2'd1, 4'd0, 3'd1, 8.5280e-3, 0.0260);  // coastal      2 m
        run_case(2'd1, 4'd2, 3'd1, 1.9413e-3, 0.0437);  // coastal      6 m
        run_case(2'd1, 4'd4, 3'd1, 4.4854e-4, 0.0596);  // coastal     10 m
        run_case(2'd2, 4'd0, 3'd1, 5.1234e-3, 0.0395);  // harbor     0.5 m
        run_case(2'd2, 4'd2, 3'd1, 4.2804e-4, 0.0627);  // harbor     1.5 m

        // ---------- TEST 1b: knob buoc song ----------
        $display("");
        $display("[TEST 1b] Knob buoc song lam DOI ket qua? (0=450 1=532 2=650 nm)");
        $display("  Ky vong theo Kebapci'23 §I: clear ocean uu blue (450),");
        $display("  nuoc duc dich dan ve red (650).");
        begin : lam_test
            integer c450, c532, c650;
            // clear ocean 15 m
            water_type = 2'd0; dist_idx = 4'd2; turb_level = 3'd1;
            lambda_sel = 2'd0; n_click=0; n_err=0; n_lost=0;
            repeat(20) @(negedge clk);
            for (i=0;i<N_QUBIT;i=i+1) pulse_qubit; c450 = n_click;
            lambda_sel = 2'd1; n_click=0; n_err=0; n_lost=0;
            repeat(20) @(negedge clk);
            for (i=0;i<N_QUBIT;i=i+1) pulse_qubit; c532 = n_click;
            lambda_sel = 2'd2; n_click=0; n_err=0; n_lost=0;
            repeat(20) @(negedge clk);
            for (i=0;i<N_QUBIT;i=i+1) pulse_qubit; c650 = n_click;
            $display("  clear ocean 15m : click 450=%0d 532=%0d 650=%0d", c450, c532, c650);
            if (c450 > c532 && c532 > c650)
                $display("    [OK] blue > green > red o clear ocean");
            else begin
                $display("    *** FAIL: thu tu buoc song sai o clear ocean");
                fails = fails + 1;
            end

            // harbor 0.5 m (nuoc rat duc)
            water_type = 2'd2; dist_idx = 4'd0;
            lambda_sel = 2'd0; n_click=0; n_err=0; n_lost=0;
            repeat(20) @(negedge clk);
            for (i=0;i<N_QUBIT;i=i+1) pulse_qubit; c450 = n_click;
            lambda_sel = 2'd2; n_click=0; n_err=0; n_lost=0;
            repeat(20) @(negedge clk);
            for (i=0;i<N_QUBIT;i=i+1) pulse_qubit; c650 = n_click;
            $display("  harbor      0.5m: click 450=%0d 650=%0d", c450, c650);
            if (c650 > c450)
                $display("    [OK] red > blue o harbor (nuoc duc)");
            else begin
                $display("    *** FAIL: red khong tot hon blue o harbor");
                fails = fails + 1;
            end
        end
        lambda_sel = 2'd0;

        // ---------- TEST 2 ----------
        $display("");
        $display("[TEST 2] Duong tin hieu OOK: khung [SYNC][b1=basis][b0=data][IDLE]");

        // 2a. Ep no_click: p_sig = 0 va p_noise = 0 -> luon mat photon
        mu_level = 4'd0; p_noise = 24'd0;
        water_type = 2'd0; dist_idx = 4'd0; turb_level = 3'd1;
        repeat (20) @(negedge clk);
        pulse_qubit;
        send_frame(1'b1, 1'b1);
        if (sync_hi == 0 && b1_seen == 0 && b0_seen == 0)
            $display("  [OK] no_click -> im lang ca khung (sync=%0d b1=%0d b0=%0d)",
                     sync_hi, b1_seen, b0_seen);
        else begin
            $display("  *** FAIL: no_click nhung van co tin hieu ra (sync=%0d b1=%0d b0=%0d)",
                     sync_hi, b1_seen, b0_seen);
            fails = fails + 1;
        end

        // 2b. Ep click KHONG loi: p_sig = max, epol khong dung vi ta ep bang mu
        //     -> dung dist gan nhat + mu max de p_sig ~ 1, va kiem tra khung
        //     truyen nguyen ven khi err_inject = 0.
        mu_level = 4'd15; p_noise = 24'd0;
        water_type = 2'd0; dist_idx = 4'd0;
        repeat (20) @(negedge clk);
        // lap den khi bat duoc mot qubit click ma khong loi
        // (P_sig ~ 2e-2 o day nen can vai nghin lan thu)
        for (k = 0; k < 20000; k = k + 1) begin
            pulse_qubit;
            if (click && !err_inject) k = 100000;
        end
        if (k >= 100000) begin
            send_frame(1'b1, 1'b1);
            if (sync_hi == slot_width && b1_seen == slot_width && b0_seen == slot_width)
                $display("  [OK] click khong loi -> khung truyen nguyen ven");
            else begin
                $display("  *** FAIL: khung bi sai (sync=%0d b1=%0d b0=%0d, ky vong %0d)",
                         sync_hi, b1_seen, b0_seen, slot_width);
                fails = fails + 1;
            end
        end else
            $display("  [SKIP] khong bat duoc qubit click-khong-loi");

        // 2c. Ep click CO loi: chon cu ly xa -> epol lon; lap den khi err_inject=1
        water_type = 2'd2; dist_idx = 4'd6;   // harbor 3.5 m -> e_pol ~ 5%
        mu_level = 4'd15; p_noise = 24'h7F_FFFF;   // ep noise_det ~50% de co click
        repeat (20) @(negedge clk);
        for (k = 0; k < 5000; k = k + 1) begin
            pulse_qubit;
            if (click && err_inject) k = 100000;
        end
        if (k >= 100000) begin
            send_frame(1'b1, 1'b1);
            // SYNC va b1 phai NGUYEN VEN, chi b0 bi dao (1 -> 0)
            if (sync_hi == slot_width && b1_seen == slot_width && b0_seen == 0)
                $display("  [OK] err_inject -> CHI dao slot du lieu b0, SYNC/basis nguyen ven");
            else begin
                $display("  *** FAIL: dao sai slot (sync=%0d b1=%0d b0=%0d, ky vong %0d/%0d/0)",
                         sync_hi, b1_seen, b0_seen, slot_width, slot_width);
                fails = fails + 1;
            end
        end else
            $display("  [SKIP] khong bat duoc qubit click-co-loi");

        // ---------- TEST 3 ----------
        $display("");
        $display("[TEST 3] Thoi gian ket hop: h phai DONG BANG giua hai coh_tick");
        chan_enable = 1'b1; mu_level = 4'd8; p_noise = 24'd218;
        water_type = 2'd0; dist_idx = 4'd8; turb_level = 3'd5;
        repeat (20) @(negedge clk);
        pulse_qubit;
        prev_hs = h_s; prev_ho = h_o; changes = 0;
        // (1<<COH)/2 qubit lien tiep, moi qubit ton 3 chu ky -> van trong 1 khoi
        for (i = 0; i < (1 << COH) / 8; i = i + 1) begin
            pulse_qubit;
            if (h_s !== prev_hs || h_o !== prev_ho) changes = changes + 1;
            prev_hs = h_s; prev_ho = h_o;
        end
        if (changes <= 1)
            $display("  [OK] h_s/h_o doi %0d lan trong %0d qubit lien tiep (<=1)",
                     changes, (1 << COH) / 8);
        else begin
            $display("  *** FAIL: h doi %0d lan -> lay mau qua nhanh, sai vat ly",
                     changes);
            fails = fails + 1;
        end

        // Nhung qua nhieu khoi ket hop thi h PHAI doi
        changes = 0; prev_hs = h_s;
        for (i = 0; i < 20; i = i + 1) begin
            repeat (1 << COH) @(negedge clk);
            pulse_qubit;
            if (h_s !== prev_hs) changes = changes + 1;
            prev_hs = h_s;
        end
        if (changes >= 5)
            $display("  [OK] qua 20 khoi ket hop, h_s doi %0d lan (>=5)", changes);
        else begin
            $display("  *** FAIL: h_s chi doi %0d lan qua 20 khoi -> LFSR/ROM sai",
                     changes);
            fails = fails + 1;
        end

        // ---------- TEST 4: bypass ----------
        $display("");
        $display("[TEST 4] chan_enable = 0 -> kenh ly tuong");
        chan_enable = 1'b0;
        repeat (20) @(negedge clk);
        n_click = 0; n_lost = 0; n_err = 0;
        for (i = 0; i < 1000; i = i + 1) pulse_qubit;
        send_frame(1'b1, 1'b0);
        if (n_click == 1000 && n_lost == 0 && n_err == 0
            && sync_hi == slot_width && b1_seen == slot_width && b0_seen == 0)
            $display("  [OK] bypass: 1000/1000 click, 0 loi, khung nguyen ven");
        else begin
            $display("  *** FAIL bypass: click=%0d lost=%0d err=%0d sync=%0d b1=%0d b0=%0d",
                     n_click, n_lost, n_err, sync_hi, b1_seen, b0_seen);
            fails = fails + 1;
        end

        // ---------- TEST 5 [v12] ----------
        // The block length must be counted in QUBIT EVENTS, not clock cycles.
        // This is the regression test for the bug that made every measured
        // turbulence level look identical: with a free-running clock timer the
        // board (≈3500 qubit/s) held h constant for only ~18 qubits instead of the
        // 50 000 the model assumes, so every measurement window averaged over
        // ~10^5 independent fading samples.
        $display("");
        $display("[TEST 5] coh_sel != 0 -> h dong bang theo SO QUBIT, khong theo dong ho");
        chan_enable = 1'b1; mu_level = 4'd8; p_noise = 24'd218;
        water_type = 2'd0; dist_idx = 4'd4; turb_level = 3'd5;
        lambda_sel = 2'd0;
        coh_sel_r = 4'd1;                       // 2^(1+5) = 64 qubit / khoi
        repeat (4) @(negedge clk);
        stats_rst_r = 1'b1; @(negedge clk); stats_rst_r = 1'b0;
        repeat (20) @(negedge clk);

        // 5a. Dong ho chay rat lau ma KHONG co qubit nao -> h phai giu nguyen.
        //     stats_rst vua dat lai qub_cnt = 0 nen hai qubit duoi day la so 1 va 2,
        //     chac chan nam trong CUNG mot khoi 64 qubit.
        pulse_qubit;
        prev_hs = h_s; prev_ho = h_o;
        repeat (300000) @(negedge clk);         // 6 ms @50 MHz, gap doi tau_coh cu
        pulse_qubit;
        if (h_s === prev_hs && h_o === prev_ho)
            $display("  [OK] 300000 chu ky khong co qubit -> h khong doi (%0d/%0d)",
                     h_s, h_o);
        else begin
            $display("  *** FAIL: h doi tu (%0d/%0d) sang (%0d/%0d) khi khong co qubit nao",
                     prev_hs, prev_ho, h_s, h_o);
            fails = fails + 1;
        end

        // 5b. qub_index phai dem MOI QUBIT, ke ca qubit KHONG click.
        //     Day la loi da lam hong mot phien do 4 gio: top_module bao ve
        //     total_cnt, ma total_cnt chi tang trong FSM_PROCESS — trang thai ma
        //     qubit mat photon KHONG BAO GIO di qua. Truong bao ve khi do dem
        //     CLICK chu khong phai ATTEMPT, nen moi click deu roi vao khoi 0.
        begin : attempt_test
            integer q0, q1, n_pulse, n_click_here;
            q0 = qub_index; n_click_here = 0;
            n_pulse = 500;
            for (i = 0; i < n_pulse; i = i + 1) begin
                pulse_qubit;
                if (click) n_click_here = n_click_here + 1;
            end
            q1 = qub_index;
            $display("  qub_index: %0d -> %0d qua %0d qubit (trong do %0d click)",
                     q0, q1, n_pulse, n_click_here);
            if (q1 - q0 == n_pulse)
                $display("  [OK] qub_index dem dung SO QUBIT, khong phai so click");
            else begin
                $display("  *** FAIL: qub_index tang %0d, ky vong %0d%s",
                         q1 - q0, n_pulse,
                         (q1 - q0 == n_click_here) ? "  (dang dem CLICK!)" : "");
                fails = fails + 1;
            end
        end

        // 5c. Khoang cach giua hai lan h doi phai la BOI SO cua 64 qubit.
        begin : coh_gap_test
            integer last_change, gap, bad_gap, n_chg;
            pulse_qubit;
            prev_hs = h_s; prev_ho = h_o;
            last_change = 0; bad_gap = 0; n_chg = 0;
            for (i = 1; i <= 640; i = i + 1) begin
                pulse_qubit;
                if (h_s !== prev_hs || h_o !== prev_ho) begin
                    n_chg = n_chg + 1;
                    // Bo qua khoang dau tien: diem xuat phat khong nam tren bien khoi.
                    if (n_chg > 1) begin
                        gap = i - last_change;
                        if ((gap % 64) != 0) begin
                            bad_gap = bad_gap + 1;
                            $display("    gap sai nhip: %0d qubit", gap);
                        end
                    end
                    last_change = i;
                end
                prev_hs = h_s; prev_ho = h_o;
            end
            $display("  640 qubit, coh_sel=1 (64 qubit/khoi): h doi %0d lan, %0d lan sai nhip",
                     n_chg, bad_gap);
            // <10 vi hai khoi lien tiep co the boc trung gia tri ROM.
            if (n_chg >= 5 && bad_gap == 0)
                $display("  [OK] moi lan doi deu roi dung bien khoi 2^(coh_sel+5) qubit");
            else begin
                $display("  *** FAIL: nhip lay mau khong khop 2^(coh_sel+5) qubit");
                fails = fails + 1;
            end
        end
        coh_sel_r = 4'd0;

        // ---------- TEST 6 [rev B] ----------
        // Ablation rows 6/7 must be selected verbatim (no range offset: dist 0
        // has offset -2 and would otherwise pull the level down) and every
        // drawn h_o must come from that ROM row with the row's spread.
        $display("");
        $display("[TEST 6] Hang ablation 6/7: dung hang ROM, bo qua offset cu ly");
        begin : ablation_test
            integer row, n_draw, n_bad, k2, found, hv;
            real s1, s2, mean_h, var_h, var_rom, m_rom;
            chan_enable = 1'b1; mu_level = 4'd8; p_noise = 24'd218;
            decoy_en = 1'b0; mu_gain = 2'd0;
            water_type = 2'd0; dist_idx = 4'd0;
            coh_sel_r = 4'd1;                       // 64 qubit / khoi
            n_draw = 4096;
            for (row = 6; row <= 7; row = row + 1) begin
                turb_level = row;
                repeat (4) @(negedge clk);
                stats_rst_r = 1'b1; @(negedge clk); stats_rst_r = 1'b0;
                repeat (4) @(negedge clk);
                if (cur_level !== row) begin
                    $display("  *** FAIL: turb_level=%0d cho current_level=%0d", row, cur_level);
                    fails = fails + 1;
                end
                // hv is a 32-bit integer on purpose: the product of two 12-bit
                // ROM words would otherwise be evaluated at 12 bits and wrap.
                s1 = 0.0; s2 = 0.0;
                for (k2 = 0; k2 < 256; k2 = k2 + 1) begin
                    hv = dut.ho_rom[row*256 + k2];
                    s1 = s1 + hv;
                    s2 = s2 + hv * hv;
                end
                m_rom   = s1 / 256.0;
                var_rom = s2 / 256.0 / (m_rom * m_rom) - 1.0;
                s1 = 0.0; s2 = 0.0; n_bad = 0;
                for (i = 0; i < n_draw; i = i + 1) begin
                    repeat (64) pulse_qubit;        // mot khoi day du
                    found = 0;
                    for (k2 = 0; k2 < 256; k2 = k2 + 1)
                        if (dut.ho_rom[row*256 + k2] == h_o) found = 1;
                    if (!found) n_bad = n_bad + 1;
                    hv = h_o;
                    s1 = s1 + hv;
                    s2 = s2 + hv * hv;
                end
                mean_h = s1 / n_draw;
                var_h  = s2 / n_draw / (mean_h * mean_h) - 1.0;
                $display("  row %0d: cur_level=%0d  E[h_o]=%.1f (ROM %.1f)  sigma2=%.3f (ROM %.3f)  ngoai hang ROM=%0d",
                         row, cur_level, mean_h, m_rom, var_h, var_rom, n_bad);
                if (n_bad != 0) begin
                    $display("    *** FAIL: h_o khong thuoc hang ROM %0d", row);
                    fails = fails + 1;
                end
                if (var_h < 0.65 * var_rom || var_h > 1.35 * var_rom) begin
                    $display("    *** FAIL: phuong sai h_o lech hang ROM %0d", row);
                    fails = fails + 1;
                end
            end
            turb_level = 3'd1;
            coh_sel_r  = 4'd0;
        end

        // ---------- TEST 7 [rev B] ----------
        // Vacuum class must never click without background, the decoy class
        // must click at signal/2^shift, and mu_gain_log2 must scale the signal.
        $display("");
        $display("[TEST 7] Lop cuong do: vacuum = 0, decoy = signal >> shift, gain = 2^g");
        begin : decoy_test
            integer n_s, n_d, n_v, n_g0, n_g2;
            real r_sd, r_g;
            chan_enable = 1'b1; water_type = 2'd0; dist_idx = 4'd0;   // 5 m
            turb_level = 3'd0;                      // no fading: clean ratios
            coh_sel_r = 4'd0; p_noise = 24'd0; mu_level = 4'd15; mu_gain = 2'd0;
            decoy_en = 1'b1; decoy_shift = 3'd3;
            repeat (20) @(negedge clk);
            att_class = 2'd0; n_click = 0;
            for (i = 0; i < 100000; i = i + 1) pulse_qubit;
            n_s = n_click;
            att_class = 2'd1; n_click = 0;
            for (i = 0; i < 100000; i = i + 1) pulse_qubit;
            n_d = n_click;
            att_class = 2'd2; n_click = 0;
            for (i = 0; i < 100000; i = i + 1) pulse_qubit;
            n_v = n_click;
            decoy_en = 1'b0; att_class = 2'd0;
            mu_level = 4'd8; mu_gain = 2'd0;
            repeat (4) @(negedge clk);
            n_click = 0;
            for (i = 0; i < 100000; i = i + 1) pulse_qubit;
            n_g0 = n_click;
            mu_gain = 2'd2;
            repeat (4) @(negedge clk);
            n_click = 0;
            for (i = 0; i < 100000; i = i + 1) pulse_qubit;
            n_g2 = n_click;
            mu_gain = 2'd0;
            r_sd = (n_d > 0) ? (n_s * 1.0 / n_d) : 0.0;
            r_g  = (n_g0 > 0) ? (n_g2 * 1.0 / n_g0) : 0.0;
            $display("  signal=%0d decoy(>>3)=%0d vacuum=%0d  ti so=%.2f (ky vong 8)", n_s, n_d, n_v, r_sd);
            $display("  gain 1: %0d  gain 4: %0d  ti so=%.2f (ky vong 4)", n_g0, n_g2, r_g);
            if (n_v != 0) begin
                $display("    *** FAIL: lop vacuum co click khi Y0 = 0");
                fails = fails + 1;
            end
            if (r_sd < 6.5 || r_sd > 9.8) begin
                $display("    *** FAIL: ti so signal/decoy sai");
                fails = fails + 1;
            end
            if (r_g < 3.6 || r_g > 4.4) begin
                $display("    *** FAIL: mu_gain_log2 khong nhan cuong do dung 2^g");
                fails = fails + 1;
            end
            p_noise = 24'd218; turb_level = 3'd1;
        end

        // ---------- TEST 8 [rev B] ----------
        // Fresh random words: with P = 1/256 from the noise comparator alone,
        // consecutive clicks must be as rare as for independent attempts. The
        // old 1-bit-per-clock LFSR windows shared most of their bits between
        // consecutive attempts and produced many times more pairs.
        $display("");
        $display("[TEST 8] Tu ngau nhien moi cho moi attempt");
        begin : fresh_test
            integer n_pairs, n_cl, prev_click;
            real exp_pairs;
            chan_enable = 1'b1; decoy_en = 1'b0; mu_gain = 2'd0;
            mu_level = 4'd0;                        // no signal
            p_noise = 24'h010000;                   // P = 1/256
            turb_level = 3'd1; water_type = 2'd0; dist_idx = 4'd4;
            repeat (20) @(negedge clk);
            n_pairs = 0; n_cl = 0; prev_click = 0;
            for (i = 0; i < 400000; i = i + 1) begin
                pulse_qubit;
                if (click) begin
                    n_cl = n_cl + 1;
                    if (prev_click) n_pairs = n_pairs + 1;
                end
                prev_click = click;
            end
            exp_pairs = n_cl * (n_cl / 400000.0);
            $display("  %0d click / 400000 (ky vong 1562), %0d cap lien tiep (doc lap ~%.1f)",
                     n_cl, n_pairs, exp_pairs);
            if (n_cl < 1400 || n_cl > 1730) begin
                $display("    *** FAIL: tan suat click noise sai");
                fails = fails + 1;
            end
            if (n_pairs > 20) begin
                $display("    *** FAIL: click lien tiep tuong quan -> tu ngau nhien khong moi");
                fails = fails + 1;
            end
            p_noise = 24'd218; mu_level = 4'd8;
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
