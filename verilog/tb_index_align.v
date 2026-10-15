// ============================================================================
// tb_index_align.v  -- verification of the [Q2] pipeline fix in uwoc_channel.v
//
// The PATCH tells you to move the whole always block at 384-418 onto
// sample_pulse_d2. That block also contains `qub_index <= qub_cnt`, and qub_cnt
// has ALREADY been incremented by the coherence block at the sample_pulse edge,
// so evaluating it two clocks later latches the NEXT attempt's index.
//
// That failure is invisible to every host-side gate: h_f_out would be delayed by
// the same two clocks, so h_min == h_max still passes and index_step_bad still
// sees a clean +1 step. Only a probe of qub_index against a ground-truth attempt
// counter catches it. This is that probe.
//
// CHECK 1  qub_index at the event == the 0-based attempt number
// CHECK 2  h_f_out is constant across each coherence block of 2^(coh_sel+5)
//          attempts and the block boundary lands exactly on the index boundary
// ============================================================================

`timescale 1ns/1ps

module tb_index_align;

    localparam COH_SEL   = 4'd1;              // 2^(1+5) = 64 attempts per block
    localparam BLK       = 64;
    localparam N_ATTEMPT = 8 * BLK;

    reg clk = 1'b0;
    reg rst_n = 1'b0;
    always #10 clk = ~clk;                    // 50 MHz

    reg  sample_en = 1'b0;
    wire click, no_click, err_inject;
    wire [11:0] h_s_out, h_o_out, h_f_out;
    wire [23:0] qub_index;
    wire [2:0]  current_level;
    wire [15:0] click_count, lost_count, flip_count;
    wire [15:0] nexp_inv;

    uwoc_channel #(.COH_LOG2(18)) dut (
        .clk(clk), .rst_n(rst_n),
        .signal_in(1'b0), .signal_out(),
        .slot_width(24'd500),
        .chan_enable(1'b1),
        .water_type(2'd0), .dist_idx(4'd4), .turb_level(3'd5),
        .lambda_sel(2'd0), .mu_level(4'd8),
        .p_noise(24'd218), .dynamic_enable(1'b0),
        .mu_gain_log2(2'd0), .decoy_enable(1'b0),
        .att_class(2'd0), .decoy_shift(3'd0),
        .coh_sel(COH_SEL), .stats_rst(1'b0),
        .current_level(current_level),
        .click(click), .no_click(no_click), .err_inject(err_inject),
        .click_count(click_count), .lost_count(lost_count),
        .flip_count(flip_count),
        .h_s_out(h_s_out), .h_o_out(h_o_out), .h_f_out(h_f_out),
        .qub_index(qub_index), .nexp_inv(nexp_inv),
        .sample_en(sample_en)
    );

    integer n_issued;         // ground truth: attempts strobed so far
    integer n_events;
    integer fails;
    integer i;
    integer idx_err, h_err;
    reg [11:0] blk_h;
    integer    blk_of_h;

    // Observe the event. `click`/`no_click` pulse two clocks after sample_pulse
    // now, so qub_index has long settled; it must still read the index of the
    // attempt that produced this event.
    always @(posedge clk) begin
        if (rst_n && (click || no_click)) begin
            n_events = n_events + 1;
            if (qub_index !== (n_events - 1)) begin
                if (idx_err < 5)
                    $display("FAIL: event %0d reports qub_index=%0d, expected %0d",
                             n_events, qub_index, n_events - 1);
                idx_err = idx_err + 1;
            end
            // h must be constant inside a coherence block of BLK attempts
            if (((n_events - 1) % BLK) == 0) begin
                blk_h    = h_f_out;
                blk_of_h = (n_events - 1) / BLK;
            end
            else if (h_f_out !== blk_h) begin
                if (h_err < 5)
                    $display("FAIL: h changed inside block %0d at index %0d (%0d -> %0d)",
                             blk_of_h, n_events - 1, blk_h, h_f_out);
                h_err = h_err + 1;
            end
        end
    end

    // Strictly serial pacing: one attempt every 5 clocks, exactly what
    // realtime_attempt_engine issues at rate_sel = 6.
    initial begin
        n_issued = 0; n_events = 0; fails = 0; idx_err = 0; h_err = 0;
        blk_h = 12'd0; blk_of_h = 0;
        repeat (6) @(posedge clk);
        rst_n = 1'b1;
        repeat (6) @(posedge clk);

        for (i = 0; i < N_ATTEMPT; i = i + 1) begin
            @(posedge clk) sample_en <= 1'b1;
            @(posedge clk) sample_en <= 1'b0;
            n_issued = n_issued + 1;
            repeat (3) @(posedge clk);
        end
        repeat (40) @(posedge clk);

        $display("");
        $display("attempts strobed : %0d", n_issued);
        $display("events observed  : %0d", n_events);
        $display("index mismatches : %0d", idx_err);
        $display("h-in-block errors: %0d", h_err);
        $display("last qub_index   : %0d", qub_index);

        if (n_events !== n_issued) begin
            $display("FAIL: %0d events for %0d attempts", n_events, n_issued);
            fails = fails + 1;
        end
        if (idx_err != 0) fails = fails + 1;
        if (h_err   != 0) fails = fails + 1;
        // A dead channel would trivially pass the index check, so make sure the
        // fading actually moved during the run.
        if (click_count == 16'd0 && lost_count == 16'd0) begin
            $display("FAIL: channel produced no decisions at all");
            fails = fails + 1;
        end

        $display("");
        if (fails == 0) $display("=== INDEX ALIGNMENT OK ===");
        else            $display("=== %0d FAILURE(S) ===", fails);
        $finish;
    end

endmodule
