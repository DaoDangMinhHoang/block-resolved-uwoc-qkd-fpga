// ============================================================
// TOP MODULE v9.0 — BB84 QKD with PC-Interactive Demo Mode
// ============================================================
//
// UPGRADE v9 (from v8.1):
//   [1] UART RX: receive Alice's bits from the PC (bb84_fpga_demo.py)
//   [2] PC_INPUT mode (SW[9]=1): Alice data from UART instead of the TRNG
//       → Flow: PC sends bits → FPGA processes → results returned to the PC
//   [3] Per-qubit response: send each qubit's result back to the PC
//   [4] Batch mode: the PC sends a bit stream, the FPGA processes it in order
//   [5] Backward compatible: SW[9]=0 → behaves exactly like v8.1
//
// DEMO FLOW (SW[9]=1, mode PC Input):
//   PC (Alice) → UART byte → FPGA decode → OOK TX → Channel →
//   OOK RX → Bob → Error Est → UART response → PC (Bob receives)
//
// UART COMMAND PROTOCOL (PC → FPGA):
//   Byte[7]=1: Qubit command
//     Bit[2] = alice_data
//     Bit[1] = alice_basis  (0=Z/Rectilinear, 1=X/Diagonal)
//     Bit[0] = bob_basis    (0=Z, 1=X)
//   Byte[7]=0: Control command
//     0x01 = Reset statistics
//     0x02 = Request status report (existing format)
//     0x10 = Echo test (FPGA replies 0xAA)
//
// UART RESPONSE (FPGA → PC, one line per CLICK — a lost photon emits nothing):
//   @a_data,a_basis,b_basis,bob_bit,bmatch,error,irrad,total,sifted,errors*\r\n
//     a_data…error : 1 decimal digit each
//     irrad        : 3 decimal digits, h on the 128 = 1.0 scale
//     total        : 6 HEX digits — [v12] widened from 4. Attempt index of THIS
//                    click, reset by 0x01. Divide by 2^(coh_sel+5) to get the
//                    coherence block the click belongs to.
//     sifted,errors: 4 hex digits each, running totals
//     mode         : [v13] 0=AGGRESSIVE 1=MODERATE 2=CONSERVATIVE 3=PAUSE
//     mu           : [v13] 1 hex digit, intensity level actually applied (8 = nominal)
//     lam          : [v13] 0=450 1=532 2=650 nm, chosen by the controller
//   Line length 42 bytes ⇒ 3.65 ms @115 200, inside FSM_REPORT's 4.4 ms.
//   mode, mu and lam are APPENDED, so fields 0..9 keep their positions and every
//   log or parser written against [v12] still reads correctly.
//
// FSM states: IDLE → WAIT_CMD → ENCODE → TX_WAIT → RX_WAIT →
//             PROCESS → REPORT → GAP → NEXT → IDLE
// ============================================================

module top_module (
    input  wire        CLOCK_50,
    input  wire [9:0]  SW,
    input  wire [3:0]  KEY,
    output wire [9:0]  LEDR,
    output wire [7:0]  LEDG,
    output wire [6:0]  HEX0,
    output wire [6:0]  HEX1,
    output wire [6:0]  HEX2,
    output wire [6:0]  HEX3,
    output wire        UART_TXD,
    input  wire        UART_RXD,
    inout  wire [35:0] GPIO_0,
    inout  wire [35:0] GPIO_1
);

    // ============================
    // Control signals
    // ============================
    wire mode_auto       = ~SW[0];
    wire adaptive_enable = SW[1];
    wire chan_enable     = SW[4];       // [v10] 0 = ideal channel (bypass)
    wire pc_input_mode   = SW[9];       // [v9] PC input mode
    wire spy_active      = ~KEY[0];
    wire manual_send     = ~KEY[1];
    wire clear_stats     = ~KEY[2];
    wire rst_n           = KEY[3];

    // ============================
    // [v10] UWOC CHANNEL CONFIGURATION
    // ============================
    // The underwater parameter space is 3-dimensional (water type × range ×
    // turbulence) + wavelength — more than the DE1's 10 switches. So the config
    // lives in registers set over UART; the switches only supply RESET DEFAULTS.
    //
    //   UART (PC → FPGA), bytes with bit[7] = 0:
    //     0x30 | dist[3:0]                → set the range index  (0x30–0x3F)
    //     0x40 | {water[1:0], lam[1:0]}   → set water type + λ   (0x40–0x4F)
    //     0x50 | turb[2:0]                → set turbulence level (0x50–0x57),
    //                                       random walk of the level OFF
    //     0x58 | turb[2:0]                → same, random walk ON    (0x58–0x5F)
    //     0x60 | coh[3:0]                 → coherence-block selector (0x60–0x6F)
    //                                       0 = legacy clock timer,
    //                                       k = 2^(k+5) QUBIT EVENTS
    reg  [1:0]  cfg_water;
    reg  [3:0]  cfg_dist;
    reg  [2:0]  cfg_turb;
    reg  [1:0]  cfg_lambda;
    reg  [3:0]  cfg_coh;
    reg         cfg_dyn;

    // [v12] The ±1 random walk of the turbulence level used to be switched on
    // implicitly by `cfg_turb >= 2`. That makes the level an UNCONTROLLED variable
    // in exactly the experiment that sweeps it: a run nominally at L2 spent part of
    // its time at L3, and the mixture inflated the measured spread of h (measured
    // CV = 0.218 at 10 m where L2 alone predicts 0.135 and a 50/50 L2/L3 mixture
    // predicts 0.210). It is now explicit and defaults to OFF.
    wire dynamic_turb = chan_enable & cfg_dyn;

    // ============================
    // FSM state declarations
    // ============================
    // Must be declared BEFORE the UART command decoder below (which uses fsm and
    // FSM_ENCODE). v9 declared them at the end of the file → a forward reference:
    // Quartus tolerates it, ModelSim errors out and creates a shadow implicit net.
    localparam FSM_IDLE     = 4'd0;
    localparam FSM_WAIT_CMD = 4'd1;
    localparam FSM_ENCODE   = 4'd2;
    localparam FSM_TX_WAIT  = 4'd3;
    localparam FSM_RX_WAIT  = 4'd4;
    localparam FSM_PROCESS  = 4'd5;
    localparam FSM_REPORT   = 4'd6;
    localparam FSM_GAP      = 4'd7;
    localparam FSM_NEXT     = 4'd8;

    reg [3:0] fsm;

    // ============================
    // GPIO
    // ============================
    wire gpio_tx_out;
    wire gpio_rx_in;
    assign GPIO_0[0]    = gpio_tx_out;
    assign gpio_rx_in   = ~GPIO_0[1];
    assign GPIO_0[35:2] = {34{1'bz}};
    assign GPIO_1       = {36{1'bz}};

    // ============================
    // Adaptive controller outputs
    // ============================
    wire [3:0]  adapt_power;
    wire [7:0]  adapt_basis_prob;
    wire [23:0] adapt_slot_width;
    wire [7:0]  adapt_rep_gap;
    wire [1:0]  adapt_mode;
    wire        adapt_tx_allowed;
    wire [7:0]  adapt_key_rate;
    wire [1:0]  adapt_lambda;        // [v10] λ in use (may be a probe candidate)
    wire [1:0]  adapt_lam_best;      // [v10] converged λ
    wire [7:0]  adapt_stale;

    localparam [23:0] FIXED_SLOT = 24'd250_000;   // 5ms per slot
	 localparam [23:0] TURBO_SLOT = 24'd500;   // 10µs instead of 5ms

    wire [23:0] active_slot_width = pc_input_mode   ? TURBO_SLOT :
												adaptive_enable ? adapt_slot_width
																	 : FIXED_SLOT;
    wire [7:0]  active_basis_prob = adaptive_enable ? adapt_basis_prob : 8'd128;
    wire [3:0]  active_power      = adaptive_enable ? adapt_power : 4'd8;
    // [v10] λ: adaptive → the controller picks it; fixed → taken from cfg_lambda
    wire [1:0]  active_lambda     = adaptive_enable ? adapt_lambda : cfg_lambda;

    reg [23:0] slot_width_latched;

    // ============================
    // [v9] UART RX — Receive commands from PC
    // ============================
    wire [7:0] uart_rx_data;
    wire       uart_rx_valid;
    wire       uart_rx_busy;

    uart_rx #(.CLK_FREQ(50000000), .BAUD(115200)) uart_rx_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .rx_in(UART_RXD),
        .rx_data(uart_rx_data),
        .rx_valid(uart_rx_valid),
        .rx_busy(uart_rx_busy)
    );

    // [v9] PC command registers
    // ⚠ [v11] QUBIT COMMANDS MUST GO THROUGH A FIFO.
    //   v9/v10 held the command in a SINGLE register (pc_cmd_ready + 3 data bits).
    //   Two bugs followed directly:
    //   [1] LOST COMMANDS. Sending back-to-back, the PC spaces bytes 86.8 µs
    //       apart (10 bits @115200) while one qubit takes ≥220 µs (TX 50 µs +
    //       rx_timeout 160 µs + FSM_GAP 10 µs @ TURBO_SLOT = 500) — and a click
    //       adds ~4 ms more in FSM_REPORT. A byte arriving before the FSM is
    //       back in FSM_WAIT_CMD overwrites its predecessor and vanishes. Only
    //       ~40% of commands ran, so MEASURED P_click came out ~2.8× below the
    //       model (3.7e-3 vs 1.04e-2, clear ocean d = 5 m), sifted bits likewise.
    //   [2] BROKEN basis_match. error_estimation uses b_basis COMBINATIONALLY
    //       all the way to FSM_PROCESS, but a byte arriving mid-qubit changes
    //       pc_cmd_bbasis ⇒ this qubit's basis is compared with the next one's.
    //   A 64-deep FIFO + latching the command once at pop fixes both.
    localparam CMD_LOG2  = 6;
    localparam CMD_DEPTH = 1 << CMD_LOG2;      // 64 commands × 3 bits

    reg  [2:0] cmd_fifo [0:CMD_DEPTH-1];
    reg  [CMD_LOG2:0] cmd_wptr, cmd_rptr;      // 1 extra bit to tell full from empty
    reg  [15:0] cmd_drop_cnt;                  // commands dropped because the FIFO was full

    wire cmd_empty = (cmd_wptr == cmd_rptr);
    wire cmd_full  = (cmd_wptr[CMD_LOG2-1:0] == cmd_rptr[CMD_LOG2-1:0]) &&
                     (cmd_wptr[CMD_LOG2]     != cmd_rptr[CMD_LOG2]);
    // One-cycle pulse: the FSM sits in FSM_WAIT_CMD for exactly 1 tick when the FIFO is non-empty.
    wire cmd_pop   = pc_input_mode && (fsm == FSM_WAIT_CMD) && !cmd_empty;

    reg        pc_cmd_data;      // Alice data bit from PC
    reg        pc_cmd_abasis;    // Alice basis from PC
    reg        pc_cmd_bbasis;    // Bob basis from PC
    reg        pc_reset_req;     // Reset request from PC
    reg        pc_status_req;    // Status request from PC

    // [v9] Command decoder
    always @(posedge CLOCK_50 or negedge rst_n) begin
        if (!rst_n) begin
            pc_cmd_data   <= 1'b0;
            pc_cmd_abasis <= 1'b0;
            pc_cmd_bbasis <= 1'b0;
            pc_reset_req  <= 1'b0;
            pc_status_req <= 1'b0;
            cmd_wptr      <= {(CMD_LOG2+1){1'b0}};
            cmd_rptr      <= {(CMD_LOG2+1){1'b0}};
            cmd_drop_cnt  <= 16'd0;
            // [v10] Reset defaults come from the switches; UART overrides afterwards
            cfg_water     <= SW[3:2];
            cfg_dist      <= 4'd0;
            cfg_turb      <= SW[7:5];
            cfg_lambda    <= 2'd0;      // 450 nm
            cfg_coh       <= 4'd0;      // legacy clock timer until the host says otherwise
            cfg_dyn       <= 1'b0;      // level random walk OFF by default
        end else begin
            // Clear one-shot signals
            pc_reset_req  <= 1'b0;
            pc_status_req <= 1'b0;

            // ---- Pop: latch the command and HOLD it for the whole qubit, until the next pop ----
            if (cmd_pop) begin
                {pc_cmd_data, pc_cmd_abasis, pc_cmd_bbasis}
                    <= cmd_fifo[cmd_rptr[CMD_LOG2-1:0]];
                cmd_rptr <= cmd_rptr + 1'b1;
            end

            if (uart_rx_valid && pc_input_mode) begin
                if (uart_rx_data[7]) begin
                    // Qubit command: bit[7]=1 → push into the FIFO
                    if (!cmd_full) begin
                        cmd_fifo[cmd_wptr[CMD_LOG2-1:0]] <= uart_rx_data[2:0];
                        cmd_wptr <= cmd_wptr + 1'b1;
                    end else begin
                        // The PC is sending faster than the FPGA can process.
                        // LEDG[7] lights up: this batch's data is UNUSABLE (missing qubits).
                        cmd_drop_cnt <= cmd_drop_cnt + 16'd1;
                    end
                end else begin
                    // Control command: bit[7]=0
                    casez (uart_rx_data)
                        // Resetting the statistics also FLUSHES the FIFO: qubits from
                        // the previous batch must not count toward the next one.
                        // Assigned after the pop branch above, so it wins — as intended.
                        8'h01:      begin
                                        pc_reset_req <= 1'b1;
                                        cmd_wptr     <= {(CMD_LOG2+1){1'b0}};
                                        cmd_rptr     <= {(CMD_LOG2+1){1'b0}};
                                        cmd_drop_cnt <= 16'd0;
                                    end
                        8'h02:      pc_status_req <= 1'b1;
                        // [v10] UWOC channel configuration
                        8'b0011_????: cfg_dist   <= uart_rx_data[3:0];
                        8'b0100_????: begin
                                        cfg_water  <= uart_rx_data[3:2];
                                        cfg_lambda <= uart_rx_data[1:0];
                                      end
                        // bit[3] now carries the random-walk flag, so 0x50–0x57
                        // means "level, held fixed" and 0x58–0x5F "level, walking".
                        8'b0101_????: begin
                                        cfg_turb <= uart_rx_data[2:0];
                                        cfg_dyn  <= uart_rx_data[3];
                                      end
                        // [v12] coherence-block selector
                        8'b0110_????: cfg_coh    <= uart_rx_data[3:0];
                        default: ;
                    endcase
                end
            end
        end
    end

    // ============================
    // TRNG
    // ============================
    wire rand_a_data, rand_a_basis, rand_b_basis;
    wire [7:0] rand_a_byte, rand_b_byte;
    wire trng_advance;

    trng_random #(.SEED(16'h0000)) gen_alice_data (
        .clk(CLOCK_50), .rst_n(rst_n), .enable(trng_advance),
        .random_bit(rand_a_data), .random_byte()
    );
    trng_random #(.SEED(16'h0000)) gen_alice_basis (
        .clk(CLOCK_50), .rst_n(rst_n), .enable(trng_advance),
        .random_bit(rand_a_basis), .random_byte(rand_a_byte)
    );
    trng_random #(.SEED(16'h0000)) gen_bob_basis (
        .clk(CLOCK_50), .rst_n(rst_n), .enable(trng_advance),
        .random_bit(rand_b_basis), .random_byte(rand_b_byte)
    );

    wire alice_basis_biased = (rand_a_byte < active_basis_prob) ? 1'b0 : 1'b1;
    wire bob_basis_biased   = (rand_b_byte < active_basis_prob) ? 1'b0 : 1'b1;

    // [v9] MUX: data source depends on mode
    //   pc_input_mode=1 → PC provides alice data/basis/bob_basis
    //   pc_input_mode=0 → original behavior (TRNG or switches)
    wire a_data_src  = pc_input_mode ? pc_cmd_data   :
                       mode_auto     ? rand_a_data    : SW[2];
    wire a_basis_src = pc_input_mode ? pc_cmd_abasis  :
                       mode_auto     ? alice_basis_biased : SW[3];
    wire b_basis_src = pc_input_mode ? pc_cmd_bbasis  :
                       mode_auto     ? bob_basis_biased   : SW[8];

    // Wires used by rest of design
    wire a_data  = a_data_src;
    wire a_basis = a_basis_src;
    wire b_basis = b_basis_src;

    // ============================
    // Alice
    // ============================
    wire [1:0] tx_qubit;
    alice alice_inst (.a(a_data), .b(a_basis), .qubit(tx_qubit));

    // ============================
    // OOK TX
    // ============================
    wire tx_serial, tx_frame_done, tx_active;
    reg  tx_start;

    ook_tx_serializer tx_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .tx_start(tx_start), .qubit_in(tx_qubit),
        .slot_width(slot_width_latched),
        .serial_out(tx_serial),
        .frame_done(tx_frame_done),
        .tx_active(tx_active)
    );

    wire tx_powered;
    pwm_power pwm_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .power_level(active_power),
        .signal_in(tx_serial), .signal_out(tx_powered)
    );
    assign gpio_tx_out = tx_powered;

    // ============================
    // [v10] Channel: UWOC
    // ============================
    wire channel_in = pc_input_mode ? tx_serial :    // [v9] PC mode: always internal
                      SW[9]        ? gpio_rx_in :    // laser mode in non-PC
                                     tx_serial;

    wire turb_out;
    wire [2:0]  turb_cur_level;
    wire        ch_click, ch_no_click, ch_err;
    wire [15:0] ch_click_cnt, ch_lost_cnt, ch_flip_cnt;
    wire [11:0] ch_h_s, ch_h_o, ch_h_f;
    wire [15:0] ch_nexp_inv;
    wire [23:0] ch_qub_index;
    wire ch_sample_en;

    // Y₀ = (dark + background)·gate. With 60 Hz dark + 200 Hz background and a
    // 50 ns gate → 1.30e-5 → round(1.30e-5 · 2^24) = 218. (matches the LinkConfig default)
    localparam [23:0] P_NOISE_DEFAULT = 24'd218;

    uwoc_channel #(.COH_LOG2(18)) chan_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .signal_in(channel_in), .signal_out(turb_out),
        .slot_width(slot_width_latched),
        .chan_enable(chan_enable),
        .water_type(cfg_water),
        .dist_idx(cfg_dist),
        .turb_level(cfg_turb),
        .lambda_sel(active_lambda),
        .mu_level(active_power),
        .p_noise(P_NOISE_DEFAULT),
        .dynamic_enable(dynamic_turb),
        .current_level(turb_cur_level),
        .click(ch_click), .no_click(ch_no_click), .err_inject(ch_err),
        .click_count(ch_click_cnt), .lost_count(ch_lost_cnt),
        .flip_count(ch_flip_cnt),
        .h_s_out(ch_h_s), .h_o_out(ch_h_o), .h_f_out(ch_h_f),
        .qub_index(ch_qub_index),
        .nexp_inv(ch_nexp_inv),
        .sample_en(ch_sample_en),
        .coh_sel(cfg_coh),
        .stats_rst(pc_reset_req)
    );

    // The UART report still uses a 1-byte "intensity": take h_f (256 = 1.0) and
    // halve it onto the 128 = 1.0 scale of the old irradiance field → keeping the
    // existing Python scripts compatible.
    wire [7:0] irrad_combined = (ch_h_f[11:1] > 11'd255) ? 8'd255 : ch_h_f[8:1];

    // ============================
    // OOK RX
    // ============================
    wire [1:0] rx_qubit;
    wire rx_valid_pulse, rx_active_w, sig_detect;
    reg  rx_valid_flag;

    ook_rx_deserializer rx_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .serial_in(turb_out),
        .slot_width(slot_width_latched),
        .qubit_out(rx_qubit),
        .qubit_valid(rx_valid_pulse),
        .rx_active(rx_active_w),
        .signal_detect(sig_detect)
    );

    // ============================
    // Bob + Error estimation
    // ============================
    wire bob_decoded;
    bob bob_inst (
        .qubit(rx_qubit), .b_prime(b_basis),
        .spy_control(spy_active), .a1(bob_decoded)
    );

    wire match_w, basis_match_w, spy_detect_w;
    reg a_data_latch, a_basis_latch;

    error_estimation err_inst (
        .a(a_data_latch), .a1(bob_decoded),
        .b(a_basis_latch), .b_prime(b_basis),
        .match(match_w), .basis_match(basis_match_w),
        .spy_detect(spy_detect_w)
    );

    // ============================
    // Channel Monitor
    // ============================
    wire [7:0] mon_qber, mon_snr, mon_photon, mon_sifted;
    wire [7:0]  mon_loss, mon_jitter;
    wire [15:0] mon_photon_cnt;
    wire        mon_window, mon_valid;
    reg        evt_done, evt_lost, evt_bmatch, evt_derr;

    // ⚠ NEXP_LOG2 MUST match log2(--window) of python/uwoc_lut_gen.py
    //   (default 65536 = 2^16), because nexp_inv is generated against that window.
    channel_monitor #(.ATTEMPT_LOG2(16), .NEXP_LOG2(16), .MIN_SIFT(16)) ch_mon (
        .clk(CLOCK_50), .rst_n(rst_n),
        .evt_qubit_done(evt_done),
        .evt_qubit_lost(evt_lost),
        .evt_basis_match(evt_bmatch),
        .evt_data_error(evt_derr),
        .nexp_inv(ch_nexp_inv),
        .enable(1'b1), .clear(clear_stats | pc_reset_req),
        .qber(mon_qber), .snr_level(mon_snr),
        .photon_rate(mon_photon), .photon_count(mon_photon_cnt),
        .sifted_rate(mon_sifted),
        .loss_rate(mon_loss), .qber_jitter(mon_jitter),
        .window_valid(mon_valid),
        .window_pulse(mon_window)
    );

    // ============================
    // Adaptive Controller
    // ============================
    adaptive_controller adapt_ctrl (
        .clk(CLOCK_50), .rst_n(rst_n),
        .qber(mon_qber), .snr_level(mon_snr),
        .photon_rate(mon_photon), .photon_count(mon_photon_cnt),
        .qber_jitter(mon_jitter), .loss_rate(mon_loss),
        .window_valid(mon_valid), .window_pulse(mon_window),
        .adaptive_enable(adaptive_enable),
        // [v14] The host's 0x01 already clears ch_mon and the channel statistics;
        // the controller has to start each measurement point clean as well,
        // otherwise a point that ended stuck in PAUSE on a probe wavelength hands
        // that state to the next point and the rest of the session is noise.
        .cfg_rst(pc_reset_req),
        .manual_power(4'd8),
        .manual_basis_prob(8'd128),
        .manual_slot_width(FIXED_SLOT),
        .manual_lambda(cfg_lambda),
        .power_level(adapt_power),
        .basis_prob_z(adapt_basis_prob),
        .slot_width_out(adapt_slot_width),
        .rep_gap(adapt_rep_gap),
        .lambda_sel(adapt_lambda), .lambda_best(adapt_lam_best),
        .mode(adapt_mode),
        .tx_allowed(adapt_tx_allowed),
        .key_rate_est(adapt_key_rate),
        .stale_windows(adapt_stale)
    );

    // ============================
    // UART TX (shared between reporter and per-qubit response)
    // ============================
    wire [7:0] uart_data_w;
    wire       uart_start_w, uart_busy_w;

    uart_tx uart_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .tx_data(uart_data_w),
        .tx_start(uart_start_w),
        .tx_out(UART_TXD),
        .tx_busy(uart_busy_w)
    );

    // ============================
    // [v9] Per-Qubit Response Reporter
    // ============================
    // [v12] total_cnt is 24-bit. The host bins each click into a coherence block by
    // (total_qubits >> (coh_sel+5)); with a 16-bit counter that index wrapped every
    // 65 536 attempts, i.e. once per block at coh_sel = 11, and the far points
    // (4 clicks per 65 536 attempts at 45 m) gave the host no way to detect the wrap.
    // 24 bits covers 1.6·10^7 attempts — longer than any single measurement point.
    reg [23:0] total_cnt;
    reg [15:0] sifted_cnt, error_cnt, key_cnt;

    // Per-qubit result registers (latched in FSM_PROCESS)
    reg        result_ready;
    reg        res_a_data, res_a_basis, res_b_basis;
    reg        res_bob_bit, res_bmatch, res_error;
    reg [7:0]  res_irrad;
    reg [23:0] res_attempt;    // [v12] attempt index of THIS qubit, from the channel
    reg [1:0]  res_mode;       // [v13] adaptive mode in force for THIS qubit
    reg [3:0]  res_mu;         // [v13] intensity in force for THIS qubit
    reg [1:0]  res_lam;        // [v13] wavelength in force for THIS qubit

    // [v9] Dual-source UART multiplexer:
    //   - Per-qubit reporter (priority when result_ready)
    //   - Status reporter (existing uart_reporter)
    wire [7:0] rpt_uart_data;
    wire       rpt_uart_start;
    wire [7:0] stat_uart_data;
    wire       stat_uart_start;

    // Per-qubit response reporter instance
    per_qubit_reporter pq_rpt (
        .clk(CLOCK_50), .rst_n(rst_n),
        .result_ready(result_ready & pc_input_mode),
        .a_data(res_a_data), .a_basis(res_a_basis), .b_basis(res_b_basis),
        .bob_bit(res_bob_bit), .basis_match(res_bmatch), .data_error(res_error),
        .irradiance(res_irrad),
        // [v12] The reported "total" field is the ATTEMPT index, not total_cnt:
        // total_cnt only advances in FSM_PROCESS, which a lost photon never
        // reaches, so it counts clicks. The host needs attempts to locate the
        // coherence block.
        .total_qubits(res_attempt),
        .total_sifted(sifted_cnt), .total_errors(error_cnt),
        .adapt_mode(res_mode), .mu_level(res_mu), .lam_sel(res_lam),
        .uart_data(rpt_uart_data),
        .uart_start(rpt_uart_start),
        .uart_busy(uart_busy_w),
        .report_done()
    );

    // Status reporter (existing)
    uart_reporter logger_inst (
        .clk(CLOCK_50), .rst_n(rst_n),
        .qber(mon_qber), .snr_level(mon_snr),
        .photon_rate(mon_photon), .sifted_rate(mon_sifted),
        .power_level(active_power),
        .basis_prob(active_basis_prob),
        .slot_width(slot_width_latched),
        .rep_gap(adapt_rep_gap),
        .adapt_mode(adapt_mode),
        .turb_level(turb_cur_level),      // [v10] effective level (range-compensated)
        .fade_active(ch_no_click),        // [v10] "fade" ⇒ photon lost
        .irradiance(irrad_combined),
        .total_qubits(total_cnt[15:0]),   // status packet keeps its 4-hex field
        .total_sifted(sifted_cnt),
        .total_errors(error_cnt),
        .window_pulse(mon_window | pc_status_req),
        .enable(~pc_input_mode),  // [v9] disable auto-reporting in PC mode, at first dont have pc_status_req
        .uart_data(stat_uart_data),
        .uart_start(stat_uart_start),
        .uart_busy(uart_busy_w)
    );

    // UART MUX: per-qubit reporter takes priority in PC mode
    assign uart_data_w  = (pc_input_mode && rpt_uart_start) ? rpt_uart_data  : stat_uart_data;
    assign uart_start_w = (pc_input_mode && rpt_uart_start) ? rpt_uart_start : stat_uart_start;

    // ============================
    // Main FSM — BB84 Protocol (v9: extended with WAIT_CMD & REPORT)
    // The state declarations sit at the top of the module (see "FSM state declarations"),
    // because the UART command decoder references fsm/FSM_ENCODE before this point.
    // ============================
    reg [26:0] timeout_cnt, gap_cnt;
    reg [3:0]  key_shift;
    reg        last_error;

    wire [26:0] rx_timeout = {slot_width_latched[22:0], 4'b0000};
    wire [24:0] base_gap   = {slot_width_latched[23:0], 1'b0};
    wire [7:0]  eff_gap    = adaptive_enable ? adapt_rep_gap : 8'd1;
    wire [26:0] gap_total  = (eff_gap <= 8'd1) ? {2'b00, base_gap} :
                             (eff_gap <= 8'd2) ? {1'b0, base_gap, 1'b0} :
                             (eff_gap <= 8'd4) ? {base_gap, 2'b00} :
                             (eff_gap <= 8'd8) ? {base_gap[23:0], 3'b000} :
                                                  {base_gap[22:0], 4'b0000};

    wire tx_permitted = adaptive_enable ? adapt_tx_allowed : 1'b1;

    reg key1_prev;
    wire key1_pulse = manual_send & ~key1_prev;

    assign ch_sample_en = (fsm == FSM_ENCODE);

    // [v9] Report timing counter
    reg [19:0] report_wait_cnt;

    always @(posedge CLOCK_50 or negedge rst_n) begin
        if (!rst_n) begin
            fsm              <= FSM_IDLE;
            tx_start         <= 1'b0;
            a_data_latch     <= 1'b0;
            a_basis_latch    <= 1'b0;
            total_cnt        <= 24'd0;
            sifted_cnt       <= 16'd0;
            error_cnt        <= 16'd0;
            key_cnt          <= 16'd0;
            key_shift        <= 4'd0;
            timeout_cnt      <= 27'd0;
            gap_cnt          <= 27'd0;
            last_error       <= 1'b0;
            key1_prev        <= 1'b0;
            rx_valid_flag    <= 1'b0;
            slot_width_latched <= FIXED_SLOT;
            evt_done         <= 1'b0;
            evt_lost         <= 1'b0;
            evt_bmatch       <= 1'b0;
            evt_derr         <= 1'b0;
            result_ready     <= 1'b0;
            res_a_data       <= 1'b0;
            res_a_basis      <= 1'b0;
            res_b_basis      <= 1'b0;
            res_bob_bit      <= 1'b0;
            res_bmatch       <= 1'b0;
            res_error        <= 1'b0;
            res_irrad        <= 8'd128;
            res_attempt      <= 24'd0;
            res_mode         <= 2'd1;
            res_mu           <= 4'd8;
            res_lam          <= 2'd0;
            report_wait_cnt  <= 20'd0;
        end else begin
            key1_prev <= manual_send;

            evt_done   <= 1'b0;
            evt_lost   <= 1'b0;
            evt_bmatch <= 1'b0;
            evt_derr   <= 1'b0;

            if (rx_valid_pulse && (fsm == FSM_TX_WAIT || fsm == FSM_RX_WAIT))
                rx_valid_flag <= 1'b1;

            if (clear_stats || pc_reset_req) begin
                total_cnt  <= 24'd0;
                sifted_cnt <= 16'd0;
                error_cnt  <= 16'd0;
                key_cnt    <= 16'd0;
                key_shift  <= 4'd0;
                last_error <= 1'b0;
            end

            case (fsm)
                FSM_IDLE: begin
                    tx_start      <= 1'b0;
                    rx_valid_flag <= 1'b0;
                    result_ready  <= 1'b0;

                    if (pc_input_mode) begin
                        // [v9] PC mode: go to WAIT_CMD
                        fsm <= FSM_WAIT_CMD;
                    end else if (tx_permitted && (mode_auto || key1_pulse)) begin
                        // Original mode: auto or manual trigger
                        slot_width_latched <= active_slot_width;
                        fsm <= FSM_ENCODE;
                    end
                end

                FSM_WAIT_CMD: begin
                    // [v9] Wait for PC to send a qubit command
                    tx_start      <= 1'b0;
                    rx_valid_flag <= 1'b0;
                    result_ready  <= 1'b0;

                    if (!pc_input_mode) begin
                        fsm <= FSM_IDLE;  // Mode changed, go back
                    end else if (!cmd_empty) begin
                        slot_width_latched <= active_slot_width;
                        fsm <= FSM_ENCODE;
                    end
                end

                FSM_ENCODE: begin
                    a_data_latch  <= a_data;
                    a_basis_latch <= a_basis;
                    rx_valid_flag <= 1'b0;
                    tx_start      <= 1'b1;
                    fsm           <= FSM_TX_WAIT;
                end

                FSM_TX_WAIT: begin
                    tx_start <= 1'b0;
                    if (tx_frame_done) begin
                        fsm         <= FSM_RX_WAIT;
                        timeout_cnt <= 27'd0;
                    end
                end

                FSM_RX_WAIT: begin
                    timeout_cnt <= timeout_cnt + 27'd1;
                    if (rx_valid_flag)
                        fsm <= FSM_PROCESS;
                    else if (timeout_cnt >= rx_timeout) begin
                        evt_lost <= 1'b1;
                        fsm      <= FSM_GAP;
                        gap_cnt  <= 27'd0;
                    end
                end

                FSM_PROCESS: begin
                    rx_valid_flag <= 1'b0;
                    total_cnt     <= total_cnt + 24'd1;
                    evt_done      <= 1'b1;

                    // [v9] Latch per-qubit results
                    res_a_data  <= a_data_latch;
                    res_a_basis <= a_basis_latch;
                    res_b_basis <= b_basis;
                    res_bob_bit <= bob_decoded;
                    res_bmatch  <= basis_match_w;
                    res_irrad   <= irrad_combined;
                    res_attempt <= ch_qub_index;
                    res_mode    <= adapt_mode;
                    res_mu      <= active_power;
                    // [v14] WAS MISSING. res_lam was declared and wired into the
                    // reporter but never assigned anywhere, so Quartus tied it to a
                    // constant and the log's lam_idx column read 0 on every click
                    // no matter which wavelength the channel actually used. That is
                    // precisely the column that would have shown the hill-climber
                    // parked on 650 nm, so the failure was invisible for a whole
                    // session's worth of runs.
                    res_lam     <= active_lambda;

                    if (basis_match_w) begin
                        sifted_cnt <= sifted_cnt + 16'd1;
                        evt_bmatch <= 1'b1;
                        if (spy_detect_w) begin
                            error_cnt  <= error_cnt + 16'd1;
                            evt_derr   <= 1'b1;
                            last_error <= 1'b1;
                            res_error  <= 1'b1;
                        end else begin
                            key_cnt    <= key_cnt + 16'd1;
                            key_shift  <= {key_shift[2:0], bob_decoded};
                            last_error <= 1'b0;
                            res_error  <= 1'b0;
                        end
                    end else begin
                        last_error <= 1'b0;
                        res_error  <= 1'b0;
                    end

                    // [v9] Go to REPORT state in PC mode, else GAP
                    if (pc_input_mode) begin
                        result_ready    <= 1'b1;
                        report_wait_cnt <= 20'd0;
                        fsm             <= FSM_REPORT;
                    end else begin
                        fsm     <= FSM_GAP;
                        gap_cnt <= 27'd0;
                    end
                end

                FSM_REPORT: begin
                    // [v9] Wait for per-qubit response to finish sending
                    report_wait_cnt <= report_wait_cnt + 1'b1;
                    // The line must be fully out before the FSM moves on. [v13]
                    // widened it from 36 to 42 bytes = 3.65 ms @115200, leaving
                    // only 0.53 ms of the old 4.0 ms budget; 220 000 cycles
                    // (4.4 ms) restores a sane margin. Costs ~5 µs per qubit on
                    // average at P_click ~ 10⁻³ — nothing.
                    if (report_wait_cnt >= 20'd220_000) begin
                        result_ready <= 1'b0;
                        fsm          <= FSM_GAP;
                        gap_cnt      <= 27'd0;
                    end
                end

                FSM_GAP: begin
                    rx_valid_flag <= 1'b0;
                    result_ready  <= 1'b0;
                    gap_cnt       <= gap_cnt + 27'd1;

                    // [v9] In PC mode, use shorter gap
                    if (pc_input_mode) begin
                        if (gap_cnt >= {3'b000, slot_width_latched})
                            fsm <= FSM_NEXT;
                    end else begin
                        if (gap_cnt >= gap_total)
                            fsm <= FSM_NEXT;
                    end
                end

                FSM_NEXT: begin
                    rx_valid_flag <= 1'b0;
                    result_ready  <= 1'b0;
                    fsm           <= FSM_IDLE;
                end

                default: fsm <= FSM_IDLE;
            endcase
        end
    end

    assign trng_advance = (fsm == FSM_NEXT);

    // ============================
    // LEDs
    // ============================
    assign LEDR[1:0] = tx_qubit;
    assign LEDR[3:2] = rx_qubit;
    assign LEDR[4]   = tx_active;
    assign LEDR[5]   = rx_active_w;
    assign LEDR[6]   = sig_detect;
    assign LEDR[7]   = basis_match_w;
    assign LEDR[8]   = spy_active;
    assign LEDR[9]   = pc_input_mode;   // [v9] Show PC mode

    assign LEDG[1:0] = adapt_mode;
    assign LEDG[2]   = adapt_tx_allowed;
    assign LEDG[3]   = chan_enable;
    assign LEDG[4]   = ch_no_click;      // [v10] photon lost (replaces fade)
    assign LEDG[5]   = mon_valid;        // [v10] window has enough samples to be trusted
    assign LEDG[6]   = adaptive_enable;
    // [v11] Lit = the command FIFO overflowed at least once ⇒ the FPGA dropped
    //       qubits ⇒ that batch's P_click and sifted-bit count are artificially low. Lower --chunk and remeasure.
    assign LEDG[7]   = (cmd_drop_cnt != 16'd0);

    // ============================
    // 7-Segment Displays
    // ============================
    seven_seg_decoder h0 (.hex_digit(total_cnt[3:0]), .segments(HEX0));

    wire [3:0] qber_disp = (mon_qber >= 8'd100) ? 4'hF :
                           (mon_qber >= 8'd50)  ? 4'hA :
                           (mon_qber >= 8'd25)  ? 4'h6 :
                           (mon_qber >= 8'd12)  ? 4'h3 :
                           (mon_qber >= 8'd4)   ? 4'h1 : 4'h0;
    seven_seg_decoder h1 (.hex_digit(qber_disp), .segments(HEX1));
    seven_seg_decoder h2 (.hex_digit(mon_snr[7:4]), .segments(HEX2));
    seven_seg_decoder h3 (.hex_digit(fsm[3:0]), .segments(HEX3)); // [v9] Show FSM state

endmodule

// ============================================================
// PER-QUBIT REPORTER — Sends result of each qubit to PC
// ============================================================
// Format: @<a_data>,<a_basis>,<b_basis>,<bob>,<bmatch>,<err>,<irrad>,<total>,<sifted>,<errors>*\r\n
// Example: @1,0,0,1,1,0,135,0042,0021,0003*\r\n
// ============================================================
module per_qubit_reporter (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        result_ready,    // Pulse: new qubit result
    input  wire        a_data,
    input  wire        a_basis,
    input  wire        b_basis,
    input  wire        bob_bit,
    input  wire        basis_match,
    input  wire        data_error,
    input  wire [7:0]  irradiance,
    input  wire [23:0] total_qubits,   // [v12] 6 hex digits — see total_cnt
    input  wire [15:0] total_sifted,
    input  wire [15:0] total_errors,
    input  wire [1:0]  adapt_mode,     // [v13] 0=AGG 1=MOD 2=CON 3=PAUSE
    input  wire [3:0]  mu_level,       // [v13] intensity actually applied
    input  wire [1:0]  lam_sel,        // [v13] wavelength actually applied
    output reg  [7:0]  uart_data,
    output reg         uart_start,
    input  wire        uart_busy,
    output reg         report_done
);

    reg [7:0]  msg_buf [0:49];
    reg [5:0]  msg_len;
    reg [5:0]  send_idx;

    localparam S_IDLE  = 2'd0;
    localparam S_BUILD = 2'd1;
    localparam S_SEND  = 2'd2;
    localparam S_WAIT  = 2'd3;

    reg [1:0] state;
    reg       result_latched;

    // Latch result on rising edge of result_ready
    reg prev_ready;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state         <= S_IDLE;
            uart_start    <= 1'b0;
            uart_data     <= 8'd0;
            send_idx      <= 6'd0;
            msg_len       <= 6'd0;
            report_done   <= 1'b0;
            prev_ready    <= 1'b0;
        end else begin
            uart_start  <= 1'b0;
            report_done <= 1'b0;
            prev_ready  <= result_ready;

            case (state)
                S_IDLE: begin
                    if (result_ready && !prev_ready)
                        state <= S_BUILD;
                end

                S_BUILD: begin
                    // Build response message
                    msg_buf[0]  <= "@";
                    msg_buf[1]  <= a_data ? "1" : "0";
                    msg_buf[2]  <= ",";
                    msg_buf[3]  <= a_basis ? "1" : "0";
                    msg_buf[4]  <= ",";
                    msg_buf[5]  <= b_basis ? "1" : "0";
                    msg_buf[6]  <= ",";
                    msg_buf[7]  <= bob_bit ? "1" : "0";
                    msg_buf[8]  <= ",";
                    msg_buf[9]  <= basis_match ? "1" : "0";
                    msg_buf[10] <= ",";
                    msg_buf[11] <= data_error ? "1" : "0";
                    msg_buf[12] <= ",";
                    // Irradiance (3 decimal digits)
                    msg_buf[13] <= dec_digit(irradiance / 100);
                    msg_buf[14] <= dec_digit((irradiance / 10) % 10);
                    msg_buf[15] <= dec_digit(irradiance % 10);
                    msg_buf[16] <= ",";
                    // Total qubits (6 hex) — [v12] widened from 4. This field is the
                    // ATTEMPT INDEX of this click; the host needs it un-wrapped to
                    // assign the click to a coherence block. 2 extra bytes per line
                    // cost ~6% of the report time and remove a silent failure mode.
                    msg_buf[17] <= hex_char(total_qubits[23:20]);
                    msg_buf[18] <= hex_char(total_qubits[19:16]);
                    msg_buf[19] <= hex_char(total_qubits[15:12]);
                    msg_buf[20] <= hex_char(total_qubits[11:8]);
                    msg_buf[21] <= hex_char(total_qubits[7:4]);
                    msg_buf[22] <= hex_char(total_qubits[3:0]);
                    msg_buf[23] <= ",";
                    // Total sifted (4 hex)
                    msg_buf[24] <= hex_char(total_sifted[15:12]);
                    msg_buf[25] <= hex_char(total_sifted[11:8]);
                    msg_buf[26] <= hex_char(total_sifted[7:4]);
                    msg_buf[27] <= hex_char(total_sifted[3:0]);
                    msg_buf[28] <= ",";
                    // Total errors (4 hex)
                    msg_buf[29] <= hex_char(total_errors[15:12]);
                    msg_buf[30] <= hex_char(total_errors[11:8]);
                    msg_buf[31] <= hex_char(total_errors[7:4]);
                    msg_buf[32] <= hex_char(total_errors[3:0]);
                    // [v13] Adaptive mode + the intensity actually applied. APPENDED
                    // at the end so fields 0..9 keep their positions and every log
                    // and parser written against [v12] still reads correctly.
                    // Without these two the adaptive-vs-fixed comparison is a black
                    // box: P_click moves, but nothing says whether the controller
                    // was in MODERATE (mu 9/15) or CONSERVATIVE (12/15), which is
                    // the entire content of the comparison.
                    msg_buf[33] <= ",";
                    msg_buf[34] <= dec_digit({6'd0, adapt_mode});
                    msg_buf[35] <= ",";
                    msg_buf[36] <= hex_char(mu_level);
                    // λ is chosen by the hill-climber when adaptive is on, so the
                    // value the HOST asked for over 0x40 is not what the channel
                    // used. Logging it is what lets the model comparison pick the
                    // right c(λ) instead of silently assuming 450 nm.
                    msg_buf[37] <= ",";
                    msg_buf[38] <= dec_digit({6'd0, lam_sel});
                    // Terminator
                    msg_buf[39] <= "*";
                    msg_buf[40] <= "\r";
                    msg_buf[41] <= "\n";

                    msg_len  <= 6'd42;
                    send_idx <= 6'd0;
                    state    <= S_SEND;
                end

                S_SEND: begin
                    if (send_idx >= msg_len) begin
                        report_done <= 1'b1;
                        state       <= S_IDLE;
                    end else if (!uart_busy) begin
                        uart_data  <= msg_buf[send_idx];
                        uart_start <= 1'b1;
                        send_idx   <= send_idx + 1'b1;
                        state      <= S_WAIT;
                    end
                end

                S_WAIT: begin
                    if (uart_busy)
                        state <= S_SEND;
                end
            endcase
        end
    end

    function [7:0] dec_digit;
        input [7:0] val;
    begin
        if (val <= 8'd9) dec_digit = "0" + val;
        else             dec_digit = "0";
    end
    endfunction

    function [7:0] hex_char;
        input [3:0] nib;
    begin
        if (nib <= 4'd9) hex_char = "0" + {4'd0, nib};
        else             hex_char = "A" + {4'd0, nib - 4'd10};
    end
    endfunction

endmodule

// ============================================================
module seven_seg_decoder (
    input  wire [3:0] hex_digit,
    output reg  [6:0] segments
);
    always @(*) begin
        case (hex_digit)
            4'h0: segments = 7'b1000000; 4'h1: segments = 7'b1111001;
            4'h2: segments = 7'b0100100; 4'h3: segments = 7'b0110000;
            4'h4: segments = 7'b0011001; 4'h5: segments = 7'b0010010;
            4'h6: segments = 7'b0000010; 4'h7: segments = 7'b1111000;
            4'h8: segments = 7'b0000000; 4'h9: segments = 7'b0010000;
            4'hA: segments = 7'b0001000; 4'hB: segments = 7'b0000011;
            4'hC: segments = 7'b1000110; 4'hD: segments = 7'b0100001;
            4'hE: segments = 7'b0000110; 4'hF: segments = 7'b0001110;
            default: segments = 7'b1111111;
        endcase
    end
endmodule
