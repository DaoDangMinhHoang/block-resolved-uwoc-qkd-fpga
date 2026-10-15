// ============================================================================
// tb_security_chain.v   (rev E)
// Self-checking testbench for security_block_monitor -> block_reporter_bin.
//
// It exists because the bug it checks for is invisible by inspection: an
// earlier revision cleared the monitor's internal `synced` flag on the same
// edge that raised block_pulse, so the reporter -- which samples one clock
// later -- latched zero into every packet and the host rejected the entire
// capture. Only a simulation catches that.
//
// [rev E] Packet v6 (144 bytes) carries per-class counters for decoy and
// vacuum attempts. The stimulus tags every attempt with a random intensity
// class and the checker rebuilds the class-1 and class-2 counters of every
// complete block from the stimulus itself, so a mis-packed field or a class
// leaking into the wrong counter fails here instead of in a 30-minute capture.
//
// Uses ATTEMPT_LOG2 = 6 (64 attempts per block) so the run is short. The logic
// under test does not depend on the parameter value.
//
//   vlog +incdir+. tb_security_chain.v security_block_monitor.v block_reporter_bin.v
//   vsim -c -do "run -all; quit -f" tb_security_chain
// ============================================================================

`timescale 1ns/1ps

module tb_security_chain;

    localparam AL2  = 6;
    localparam WIN  = (1 << AL2);
    localparam PKT  = 144;
    localparam NBLK = 8;                 // blocks tracked by the checker

    reg clk = 1'b0;
    reg rst_n = 1'b0;
    always #5 clk = ~clk;

    // ---- stimulus ---------------------------------------------------------
    reg         ev_done, ev_lost, ev_bmatch, ev_bz, ev_derr;
    reg  [1:0]  ev_class;
    reg  [11:0] h_f;
    reg  [23:0] qub_index;
    reg  [79:0] cfg_bus;

    // ---- monitor ----------------------------------------------------------
    wire        blk_pulse, blk_synced, cfg_chg;
    wire [31:0] blk_id, perr;
    wire [17:0] att, cl_all, lo_all, tz, cz, ez, tx, cx, ex;
    wire [143:0] dcnt, vcnt;
    wire [11:0] hmean, hmin, hmax;
    wire [79:0] blkcfg;

    security_block_monitor #(.ATTEMPT_LOG2(AL2)) dut_mon (
        .clk(clk), .rst_n(rst_n), .enable(1'b1), .clear(1'b0),
        .evt_qubit_done(ev_done), .evt_qubit_lost(ev_lost),
        .evt_basis_match(ev_bmatch), .evt_basis_is_z(ev_bz),
        .evt_data_error(ev_derr), .evt_class(ev_class),
        .h_f_sample(h_f), .qub_index(qub_index), .cfg_bus(cfg_bus),
        .block_pulse(blk_pulse), .block_synced(blk_synced), .block_id(blk_id),
        .attempts(att), .clicks_all(cl_all), .lost_all(lo_all),
        .trials_z(tz), .clicks_z(cz), .errors_z(ez),
        .trials_x(tx), .clicks_x(cx), .errors_x(ex),
        .decoy_counts(dcnt), .vacuum_counts(vcnt),
        .h_mean(hmean), .h_min(hmin), .h_max(hmax),
        .blk_cfg(blkcfg), .config_changed(cfg_chg), .protocol_errors(perr)
    );

    // ---- reporter + a uart_tx that mimics the real registered-busy timing --
    wire [7:0] u_data;
    wire       u_start;
    reg        u_busy;
    wire       rpt_busy;
    wire [31:0] dropped, wdaborts;

    // Keep this equal to AL2. It is written as a sized literal rather than
    // AL2[7:0]: AL2 is an untyped integer localparam, and a part-select on one
    // is handled inconsistently across Quartus 13.0, ModelSim ASE and iverilog.
    block_reporter_bin #(.ATTEMPT_LOG2_BYTE(8'd6)) dut_rpt (
        .clk(clk), .rst_n(rst_n), .clear(1'b0),
        .block_pulse(blk_pulse), .block_synced(blk_synced),
        .security_fast_mode(1'b1), .adaptive_enable(1'b0),
        .realtime_paced(1'b1), .rate_valid(1'b1), .config_changed(cfg_chg),
        .block_id(blk_id), .attempts(att),
        .clicks_all(cl_all), .lost_all(lo_all),
        .trials_z(tz), .clicks_z(cz), .errors_z(ez),
        .trials_x(tx), .clicks_x(cx), .errors_x(ex),
        .decoy_counts(dcnt), .vacuum_counts(vcnt),
        .h_mean(hmean), .h_min(hmin), .h_max(hmax),
        .blk_cfg(blkcfg), .attempts_per_sec(32'd10000000),
        .rate_hz(32'd10000000), .protocol_errors(perr),
        .engine_timeouts(32'd0), .engine_overruns(32'd0),
        .uart_data(u_data), .uart_start(u_start), .uart_busy(u_busy),
        .reporter_busy(rpt_busy), .dropped_blocks(dropped),
        .watchdog_aborts(wdaborts)
    );

    // uart_tx model: busy is REGISTERED, so it rises one clock after start.
    // This is the exact behaviour that broke the naive handshake.
    integer bitcnt;
    reg [7:0] captured [0:4095];
    integer ncap;

    initial begin
        u_busy = 1'b0;
        bitcnt = 0;
        ncap = 0;
    end

    always @(posedge clk) begin
        if (u_start && !u_busy) begin
            u_busy  <= 1'b1;
            bitcnt  <= 2;              // shortened bit time
            captured[ncap] <= u_data;
            ncap <= ncap + 1;
        end
        else if (u_busy) begin
            if (bitcnt > 0) bitcnt <= bitcnt - 1;
            else u_busy <= 1'b0;
        end
    end

    // ---- expected per-class counters, rebuilt from the stimulus ------------
    // flat index: (class-1)*NBLK*8 + block*8 + field
    // fields: 0 attempts 1 clicks_all 2 trials_z 3 clicks_z 4 errors_z
    //         5 trials_x 6 clicks_x 7 errors_x
    integer exp_cnt [0:2*NBLK*8-1];

    task bump(input integer cc, input integer bb, input integer ff);
        begin
            if (bb < NBLK)
                exp_cnt[(cc-1)*NBLK*8 + bb*8 + ff] =
                    exp_cnt[(cc-1)*NBLK*8 + bb*8 + ff] + 1;
        end
    endtask

    // ---- drive attempts ---------------------------------------------------
    integer i;
    integer seed;
    integer blk;
    reg is_click, is_match, is_z, is_err;
    reg [1:0] cls;

    task do_attempt(input [23:0] idx, input [11:0] hval);
        begin
            seed = seed * 1103515245 + 12345;
            is_click = ((seed >> 16) % 100) < 30;    // 30% click, exaggerated
            seed = seed * 1103515245 + 12345;
            is_match = ((seed >> 16) % 2) == 0;
            seed = seed * 1103515245 + 12345;
            is_z = ((seed >> 16) % 2) == 0;
            seed = seed * 1103515245 + 12345;
            is_err = ((seed >> 16) % 10) == 0;
            seed = seed * 1103515245 + 12345;
            cls = (seed >> 16) % 3;

            blk = idx / WIN;
            if (cls != 2'd0) begin
                bump(cls, blk, 0);
                if (is_click) bump(cls, blk, 1);
                if (is_match && is_z) begin
                    bump(cls, blk, 2);
                    if (is_click) begin
                        bump(cls, blk, 3);
                        if (is_err) bump(cls, blk, 4);
                    end
                end
                if (is_match && !is_z) begin
                    bump(cls, blk, 5);
                    if (is_click) begin
                        bump(cls, blk, 6);
                        if (is_err) bump(cls, blk, 7);
                    end
                end
            end

            @(negedge clk);
            qub_index <= idx;
            h_f       <= hval;
            ev_done   <= is_click;
            ev_lost   <= ~is_click;
            ev_bmatch <= is_match;
            ev_bz     <= is_z;
            ev_derr   <= is_click & is_match & is_err;
            ev_class  <= cls;
            @(negedge clk);
            ev_done <= 1'b0; ev_lost <= 1'b0;
            ev_bmatch <= 1'b0; ev_bz <= 1'b0; ev_derr <= 1'b0;
            ev_class <= 2'd0;
            // idle clocks so the reporter has time between attempts
            repeat (14) @(negedge clk);
        end
    endtask

    // ---- checks -----------------------------------------------------------
    integer blocks_seen;
    integer fails;
    reg [7:0] chk;
    reg [7:0] flags_b;
    reg [31:0] att_w;
    integer j, p, c, f, base, val, expv;

    initial begin
        seed = 32'hC0FFEE;
        ev_done = 0; ev_lost = 0; ev_bmatch = 0; ev_bz = 0; ev_derr = 0;
        ev_class = 2'd0;
        h_f = 12'd256; qub_index = 24'd0;
        cfg_bus = 80'h0;
        cfg_bus[1:0]   = 2'd0;      // water
        cfg_bus[5:2]   = 4'd4;      // dist
        cfg_bus[8:6]   = 3'd5;      // turb
        cfg_bus[44:41] = 4'd11;     // coh_sel
        cfg_bus[72]    = 1'b1;      // decoy enable
        cfg_bus[74:73] = 2'd2;      // shift code
        cfg_bus[76:75] = 2'd1;      // class preset
        cfg_bus[78:77] = 2'd3;      // mu gain
        blocks_seen = 0;
        fails = 0;
        for (i = 0; i < 2*NBLK*8; i = i + 1) exp_cnt[i] = 0;

        repeat (4) @(negedge clk);
        rst_n = 1'b1;
        repeat (4) @(negedge clk);

        // Start mid-block on purpose: the monitor must refuse to report a
        // partial block and only arm at the next window boundary.
        for (i = WIN/2; i < 6*WIN; i = i + 1)
            do_attempt(i[23:0], 12'd200 + ((i / WIN) * 12'd40));

        repeat (2000) @(negedge clk);

        $display("");
        $display("captured %0d bytes, %0d complete packets",
                 ncap, ncap / PKT);
        if (ncap % PKT != 0) begin
            $display("FAIL: byte count %0d is not a multiple of %0d", ncap, PKT);
            fails = fails + 1;
        end

        for (p = 0; p < ncap / PKT; p = p + 1) begin
            base = p * PKT;
            chk = 8'd0;
            for (j = 0; j < PKT-1; j = j + 1)
                chk = chk ^ captured[base + j];

            flags_b = captured[base+3];
            att_w   = {captured[base+11], captured[base+10],
                       captured[base+9],  captured[base+8]};

            $display("packet %0d: magic %02X %02X ver %0d flags %02X attempts %0d decoy %02X chk %s",
                     p, captured[base+0], captured[base+1],
                     captured[base+2], flags_b, att_w, captured[base+93],
                     (chk == captured[base+PKT-1]) ? "ok" : "BAD");

            if (captured[base+0] !== 8'hA5 || captured[base+1] !== 8'h5A) begin
                $display("FAIL: bad magic in packet %0d", p);
                fails = fails + 1;
            end
            if (captured[base+2] !== 8'd6) begin
                $display("FAIL: bad version in packet %0d", p);
                fails = fails + 1;
            end
            if (chk !== captured[base+PKT-1]) begin
                $display("FAIL: checksum mismatch in packet %0d", p);
                fails = fails + 1;
            end
            // flags bit 4 is block_synced. This is the regression that matters.
            if (flags_b[4] !== 1'b1) begin
                $display("FAIL: block_synced = 0 in packet %0d -- reporter latched a stale internal flag", p);
                fails = fails + 1;
            end
            if (flags_b[7] !== 1'b1) begin
                $display("FAIL: decoy_enable flag missing in packet %0d", p);
                fails = fails + 1;
            end
            if (captured[base+93] !== cfg_bus[79:72]) begin
                $display("FAIL: decoy byte %02X, expected %02X in packet %0d",
                         captured[base+93], cfg_bus[79:72], p);
                fails = fails + 1;
            end
            if (att_w !== WIN) begin
                $display("FAIL: attempts != %0d in packet %0d", WIN, p);
                fails = fails + 1;
            end

            // Packet p reports block p+1: the first half block is never reported.
            for (c = 1; c <= 2; c = c + 1)
                for (f = 0; f < 8; f = f + 1) begin
                    j    = base + 95 + (c-1)*24 + f*3;
                    val  = captured[j] | (captured[j+1] << 8) | (captured[j+2] << 16);
                    expv = exp_cnt[(c-1)*NBLK*8 + (p+1)*8 + f];
                    if (val != expv) begin
                        $display("FAIL: packet %0d class %0d field %0d = %0d, expected %0d",
                                 p, c, f, val, expv);
                        fails = fails + 1;
                    end
                end
            blocks_seen = blocks_seen + 1;
        end

        if (blocks_seen < 2) begin
            $display("FAIL: expected at least 2 complete blocks, got %0d", blocks_seen);
            fails = fails + 1;
        end
        if (dropped !== 32'd0) begin
            $display("FAIL: dropped_blocks = %0d", dropped);
            fails = fails + 1;
        end
        if (wdaborts !== 32'd0) begin
            $display("FAIL: watchdog_aborts = %0d", wdaborts);
            fails = fails + 1;
        end
        if (perr !== 32'd0) begin
            $display("FAIL: protocol_errors = %0d", perr);
            fails = fails + 1;
        end

        $display("");
        if (fails == 0) $display("=== ALL CHECKS PASSED ===");
        else            $display("=== %0d FAILURE(S) ===", fails);
        $finish;
    end

endmodule
