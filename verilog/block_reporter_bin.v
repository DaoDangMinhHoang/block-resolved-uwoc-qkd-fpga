// ============================================================================
// block_reporter_bin.v   (rev E)
// One 144-byte little-endian binary packet per security block (packet v6).
//
// UART HANDSHAKE. uart_tx REGISTERS tx_busy, so busy rises one clock AFTER
// tx_start is sampled. A reporter that only tests (!uart_busy) fires again on
// the next clock and drops every second byte. States:
//   ST_SEND      pulse uart_start, do NOT advance idx
//   ST_WAIT_HIGH wait for uart_busy == 1
//   ST_WAIT_LOW  wait for uart_busy == 0, then advance
// The watchdog runs in ALL THREE. busy stuck high strands ST_SEND, stuck low
// strands ST_WAIT_HIGH, never falling strands ST_WAIT_LOW; one abort path
// covers all of them and increments watchdog_aborts.
//
// The watchdog only rescues the FSM. It does not save the packet. To keep the
// packet, top_module must hold UART ownership until the reporter is finished:
//     wire security_uart_owner = security_fast_mode | sec_reporter_busy;
// otherwise a 0x70 arriving mid-packet cuts the wire.
//
// PACKET v6 (144 bytes, little endian). Bytes 0..92 keep the v5 offsets.
//   0..1   A5 5A     2 version=6     3 flags
//          flags: b0 fast_mode, b1 adaptive, b2 has_sift, b3 config_changed,
//                 b4 block_synced, b5 realtime_paced, b6 rate_valid,
//                 b7 decoy_enable
//   4..7   block_id              8..11  attempts
//   12..15 clicks_all            16..19 lost_all
//   20..23 trials_z   24..27 clicks_z   28..31 errors_z
//   32..35 trials_x   36..39 clicks_x   40..43 errors_x
//   44..45 h_mean     46..47 h_min      48..49 h_max
//   50 water  51 dist_idx  52 turb  53 lambda  54 mu  55 adapt_mode
//   56..58 p_noise    59 coh_sel        60..63 slot_width
//   64..67 attempts_per_sec      68..71 rate_hz_requested
//   72..75 dropped_blocks        76..79 protocol_errors
//   80..83 engine_timeouts       84..87 engine_overruns
//   88..91 watchdog_aborts
//   92     attempt_log2
//   93     decoy byte {0, mu_gain_log2[1:0], preset[1:0], shift_code[1:0], enable}
//   94     reserved
//   95..118  class 1 (decoy)  counters, 3 bytes each, in this order:
//            attempts, clicks_all, trials_z, clicks_z, errors_z,
//            trials_x, clicks_x, errors_x
//   119..142 class 2 (vacuum) counters, same order
//   143    XOR checksum of bytes 0..142
// The totals in bytes 8..43 count every attempt whatever its class; the host
// obtains class 0 (signal) by subtraction.
//
// At 460 800 baud a v6 packet takes 144 × 21.7 µs = 3.13 ms, 48 % of the
// 6.55 ms block at 10 MHz, so no block can arrive while a packet is going out.
// ============================================================================

module block_reporter_bin #(
    parameter [7:0]  ATTEMPT_LOG2_BYTE = 8'd16,
    parameter [15:0] WATCHDOG_CYCLES   = 16'd20000
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire        clear,          // per-campaign diagnostics reset

    input  wire        block_pulse,
    input  wire        block_synced,
    input  wire        security_fast_mode,
    input  wire        adaptive_enable,
    input  wire        realtime_paced,
    input  wire        rate_valid,
    input  wire        config_changed,

    input  wire [31:0] block_id,
    input  wire [17:0] attempts,
    input  wire [17:0] clicks_all,
    input  wire [17:0] lost_all,
    input  wire [17:0] trials_z,
    input  wire [17:0] clicks_z,
    input  wire [17:0] errors_z,
    input  wire [17:0] trials_x,
    input  wire [17:0] clicks_x,
    input  wire [17:0] errors_x,
    input  wire [143:0] decoy_counts,
    input  wire [143:0] vacuum_counts,

    input  wire [11:0] h_mean,
    input  wire [11:0] h_min,
    input  wire [11:0] h_max,

    input  wire [79:0] blk_cfg,
    input  wire [31:0] attempts_per_sec,
    input  wire [31:0] rate_hz,
    input  wire [31:0] protocol_errors,
    input  wire [31:0] engine_timeouts,
    input  wire [31:0] engine_overruns,

    output reg  [7:0]  uart_data,
    output reg         uart_start,
    input  wire        uart_busy,

    output reg         reporter_busy,
    output reg  [31:0] dropped_blocks,
    output reg  [31:0] watchdog_aborts
);

    localparam ST_IDLE      = 2'd0;
    localparam ST_SEND      = 2'd1;
    localparam ST_WAIT_HIGH = 2'd2;
    localparam ST_WAIT_LOW  = 2'd3;

    localparam [7:0] LAST_BYTE = 8'd143;

    reg [1:0]  state;
    reg [7:0]  idx;
    reg [7:0]  chk;
    reg [15:0] wdog;

    reg [7:0]  r_flags;
    reg [31:0] r_block_id, r_attempts, r_clicks_all, r_lost_all;
    reg [31:0] r_tz, r_cz, r_ez, r_tx, r_cx, r_ex;
    reg [143:0] r_dec, r_vac;
    reg [15:0] r_hmean, r_hmin, r_hmax;
    reg [79:0] r_cfg;
    reg [31:0] r_rate_meas, r_rate_req, r_dropped, r_perr, r_tmo, r_ovr, r_wda;

    wire [7:0]  b_water  = {6'd0, r_cfg[1:0]};
    wire [7:0]  b_dist   = {4'd0, r_cfg[5:2]};
    wire [7:0]  b_turb   = {5'd0, r_cfg[8:6]};
    wire [7:0]  b_lambda = {6'd0, r_cfg[10:9]};
    wire [7:0]  b_mu     = {4'd0, r_cfg[14:11]};
    wire [7:0]  b_mode   = {6'd0, r_cfg[16:15]};
    wire [23:0] b_pnoise = r_cfg[40:17];
    wire [7:0]  b_coh    = {4'd0, r_cfg[44:41]};
    wire [31:0] b_slot   = {8'd0, r_cfg[68:45]};
    wire [7:0]  b_decoy  = r_cfg[79:72];

    function [7:0] pkt;
        input [7:0] k;
        begin
            case (k)
                8'd0:   pkt = 8'hA5;
                8'd1:   pkt = 8'h5A;
                8'd2:   pkt = 8'd6;
                8'd3:   pkt = r_flags;

                8'd4:   pkt = r_block_id[7:0];
                8'd5:   pkt = r_block_id[15:8];
                8'd6:   pkt = r_block_id[23:16];
                8'd7:   pkt = r_block_id[31:24];

                8'd8:   pkt = r_attempts[7:0];
                8'd9:   pkt = r_attempts[15:8];
                8'd10:  pkt = r_attempts[23:16];
                8'd11:  pkt = r_attempts[31:24];

                8'd12:  pkt = r_clicks_all[7:0];
                8'd13:  pkt = r_clicks_all[15:8];
                8'd14:  pkt = r_clicks_all[23:16];
                8'd15:  pkt = r_clicks_all[31:24];

                8'd16:  pkt = r_lost_all[7:0];
                8'd17:  pkt = r_lost_all[15:8];
                8'd18:  pkt = r_lost_all[23:16];
                8'd19:  pkt = r_lost_all[31:24];

                8'd20:  pkt = r_tz[7:0];
                8'd21:  pkt = r_tz[15:8];
                8'd22:  pkt = r_tz[23:16];
                8'd23:  pkt = r_tz[31:24];

                8'd24:  pkt = r_cz[7:0];
                8'd25:  pkt = r_cz[15:8];
                8'd26:  pkt = r_cz[23:16];
                8'd27:  pkt = r_cz[31:24];

                8'd28:  pkt = r_ez[7:0];
                8'd29:  pkt = r_ez[15:8];
                8'd30:  pkt = r_ez[23:16];
                8'd31:  pkt = r_ez[31:24];

                8'd32:  pkt = r_tx[7:0];
                8'd33:  pkt = r_tx[15:8];
                8'd34:  pkt = r_tx[23:16];
                8'd35:  pkt = r_tx[31:24];

                8'd36:  pkt = r_cx[7:0];
                8'd37:  pkt = r_cx[15:8];
                8'd38:  pkt = r_cx[23:16];
                8'd39:  pkt = r_cx[31:24];

                8'd40:  pkt = r_ex[7:0];
                8'd41:  pkt = r_ex[15:8];
                8'd42:  pkt = r_ex[23:16];
                8'd43:  pkt = r_ex[31:24];

                8'd44:  pkt = r_hmean[7:0];
                8'd45:  pkt = r_hmean[15:8];
                8'd46:  pkt = r_hmin[7:0];
                8'd47:  pkt = r_hmin[15:8];
                8'd48:  pkt = r_hmax[7:0];
                8'd49:  pkt = r_hmax[15:8];

                8'd50:  pkt = b_water;
                8'd51:  pkt = b_dist;
                8'd52:  pkt = b_turb;
                8'd53:  pkt = b_lambda;
                8'd54:  pkt = b_mu;
                8'd55:  pkt = b_mode;

                8'd56:  pkt = b_pnoise[7:0];
                8'd57:  pkt = b_pnoise[15:8];
                8'd58:  pkt = b_pnoise[23:16];
                8'd59:  pkt = b_coh;

                8'd60:  pkt = b_slot[7:0];
                8'd61:  pkt = b_slot[15:8];
                8'd62:  pkt = b_slot[23:16];
                8'd63:  pkt = b_slot[31:24];

                8'd64:  pkt = r_rate_meas[7:0];
                8'd65:  pkt = r_rate_meas[15:8];
                8'd66:  pkt = r_rate_meas[23:16];
                8'd67:  pkt = r_rate_meas[31:24];

                8'd68:  pkt = r_rate_req[7:0];
                8'd69:  pkt = r_rate_req[15:8];
                8'd70:  pkt = r_rate_req[23:16];
                8'd71:  pkt = r_rate_req[31:24];

                8'd72:  pkt = r_dropped[7:0];
                8'd73:  pkt = r_dropped[15:8];
                8'd74:  pkt = r_dropped[23:16];
                8'd75:  pkt = r_dropped[31:24];

                8'd76:  pkt = r_perr[7:0];
                8'd77:  pkt = r_perr[15:8];
                8'd78:  pkt = r_perr[23:16];
                8'd79:  pkt = r_perr[31:24];

                8'd80:  pkt = r_tmo[7:0];
                8'd81:  pkt = r_tmo[15:8];
                8'd82:  pkt = r_tmo[23:16];
                8'd83:  pkt = r_tmo[31:24];

                8'd84:  pkt = r_ovr[7:0];
                8'd85:  pkt = r_ovr[15:8];
                8'd86:  pkt = r_ovr[23:16];
                8'd87:  pkt = r_ovr[31:24];

                8'd88:  pkt = r_wda[7:0];
                8'd89:  pkt = r_wda[15:8];
                8'd90:  pkt = r_wda[23:16];
                8'd91:  pkt = r_wda[31:24];

                8'd92:  pkt = ATTEMPT_LOG2_BYTE;
                8'd93:  pkt = b_decoy;

                // ---- class 1 (decoy): 8 fields x 3 bytes ----
                8'd95:  pkt = r_dec[  0 +: 8];
                8'd96:  pkt = r_dec[  8 +: 8];
                8'd97:  pkt = {6'd0, r_dec[ 16 +: 2]};
                8'd98:  pkt = r_dec[ 18 +: 8];
                8'd99:  pkt = r_dec[ 26 +: 8];
                8'd100: pkt = {6'd0, r_dec[ 34 +: 2]};
                8'd101: pkt = r_dec[ 36 +: 8];
                8'd102: pkt = r_dec[ 44 +: 8];
                8'd103: pkt = {6'd0, r_dec[ 52 +: 2]};
                8'd104: pkt = r_dec[ 54 +: 8];
                8'd105: pkt = r_dec[ 62 +: 8];
                8'd106: pkt = {6'd0, r_dec[ 70 +: 2]};
                8'd107: pkt = r_dec[ 72 +: 8];
                8'd108: pkt = r_dec[ 80 +: 8];
                8'd109: pkt = {6'd0, r_dec[ 88 +: 2]};
                8'd110: pkt = r_dec[ 90 +: 8];
                8'd111: pkt = r_dec[ 98 +: 8];
                8'd112: pkt = {6'd0, r_dec[106 +: 2]};
                8'd113: pkt = r_dec[108 +: 8];
                8'd114: pkt = r_dec[116 +: 8];
                8'd115: pkt = {6'd0, r_dec[124 +: 2]};
                8'd116: pkt = r_dec[126 +: 8];
                8'd117: pkt = r_dec[134 +: 8];
                8'd118: pkt = {6'd0, r_dec[142 +: 2]};

                // ---- class 2 (vacuum): 8 fields x 3 bytes ----
                8'd119: pkt = r_vac[  0 +: 8];
                8'd120: pkt = r_vac[  8 +: 8];
                8'd121: pkt = {6'd0, r_vac[ 16 +: 2]};
                8'd122: pkt = r_vac[ 18 +: 8];
                8'd123: pkt = r_vac[ 26 +: 8];
                8'd124: pkt = {6'd0, r_vac[ 34 +: 2]};
                8'd125: pkt = r_vac[ 36 +: 8];
                8'd126: pkt = r_vac[ 44 +: 8];
                8'd127: pkt = {6'd0, r_vac[ 52 +: 2]};
                8'd128: pkt = r_vac[ 54 +: 8];
                8'd129: pkt = r_vac[ 62 +: 8];
                8'd130: pkt = {6'd0, r_vac[ 70 +: 2]};
                8'd131: pkt = r_vac[ 72 +: 8];
                8'd132: pkt = r_vac[ 80 +: 8];
                8'd133: pkt = {6'd0, r_vac[ 88 +: 2]};
                8'd134: pkt = r_vac[ 90 +: 8];
                8'd135: pkt = r_vac[ 98 +: 8];
                8'd136: pkt = {6'd0, r_vac[106 +: 2]};
                8'd137: pkt = r_vac[108 +: 8];
                8'd138: pkt = r_vac[116 +: 8];
                8'd139: pkt = {6'd0, r_vac[124 +: 2]};
                8'd140: pkt = r_vac[126 +: 8];
                8'd141: pkt = r_vac[134 +: 8];
                8'd142: pkt = {6'd0, r_vac[142 +: 2]};

                default: pkt = 8'd0;      // 94 reserved
            endcase
        end
    endfunction

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state           <= ST_IDLE;
            idx             <= 8'd0;
            chk             <= 8'd0;
            wdog            <= 16'd0;
            uart_data       <= 8'd0;
            uart_start      <= 1'b0;
            reporter_busy   <= 1'b0;
            dropped_blocks  <= 32'd0;
            watchdog_aborts <= 32'd0;
        end
        else if (clear) begin
            // Cumulative diagnostics are per campaign. Without this a single
            // early dropped block would make capture_valid fail for every
            // later run until the board was power cycled.
            state           <= ST_IDLE;
            idx             <= 8'd0;
            chk             <= 8'd0;
            wdog            <= 16'd0;
            uart_start      <= 1'b0;
            reporter_busy   <= 1'b0;
            dropped_blocks  <= 32'd0;
            watchdog_aborts <= 32'd0;
        end
        else begin
            uart_start <= 1'b0;

            // A block that arrives while a packet is still going out is lost.
            // block_id is transmitted, so gaps are visible host-side too.
            if (block_pulse && (state != ST_IDLE))
                dropped_blocks <= dropped_blocks + 32'd1;

            if (state == ST_IDLE) begin
                reporter_busy <= 1'b0;
                wdog          <= 16'd0;
                if (block_pulse) begin
                    r_flags <= {blk_cfg[72], rate_valid, realtime_paced, block_synced,
                                config_changed, ((clicks_z + clicks_x) != 0),
                                adaptive_enable, security_fast_mode};
                    r_block_id   <= block_id;
                    r_attempts   <= {14'd0, attempts};
                    r_clicks_all <= {14'd0, clicks_all};
                    r_lost_all   <= {14'd0, lost_all};
                    r_tz <= {14'd0, trials_z};
                    r_cz <= {14'd0, clicks_z};
                    r_ez <= {14'd0, errors_z};
                    r_tx <= {14'd0, trials_x};
                    r_cx <= {14'd0, clicks_x};
                    r_ex <= {14'd0, errors_x};
                    r_dec <= decoy_counts;
                    r_vac <= vacuum_counts;
                    r_hmean <= {4'd0, h_mean};
                    r_hmin  <= {4'd0, h_min};
                    r_hmax  <= {4'd0, h_max};
                    r_cfg       <= blk_cfg;
                    r_rate_meas <= attempts_per_sec;
                    r_rate_req  <= rate_hz;
                    r_dropped   <= dropped_blocks;
                    r_perr      <= protocol_errors;
                    r_tmo       <= engine_timeouts;
                    r_ovr       <= engine_overruns;
                    r_wda       <= watchdog_aborts;

                    idx           <= 8'd0;
                    chk           <= 8'd0;
                    reporter_busy <= 1'b1;
                    state         <= ST_SEND;
                end
            end
            else begin
                // One watchdog covering ST_SEND, ST_WAIT_HIGH and ST_WAIT_LOW.
                // It is reset on every state transition below.
                wdog <= wdog + 16'd1;

                if (wdog >= WATCHDOG_CYCLES) begin
                    idx             <= 8'd0;
                    chk             <= 8'd0;
                    wdog            <= 16'd0;
                    reporter_busy   <= 1'b0;
                    state           <= ST_IDLE;
                    watchdog_aborts <= watchdog_aborts + 32'd1;
                end
                else begin
                    case (state)
                        ST_SEND: begin
                            if (!uart_busy) begin
                                if (idx == LAST_BYTE) begin
                                    uart_data <= chk;
                                end
                                else begin
                                    uart_data <= pkt(idx);
                                    chk       <= chk ^ pkt(idx);
                                end
                                uart_start <= 1'b1;
                                wdog       <= 16'd0;
                                state      <= ST_WAIT_HIGH;
                            end
                        end

                        ST_WAIT_HIGH: begin
                            if (uart_busy) begin
                                wdog  <= 16'd0;
                                state <= ST_WAIT_LOW;
                            end
                        end

                        ST_WAIT_LOW: begin
                            if (!uart_busy) begin
                                wdog <= 16'd0;
                                if (idx == LAST_BYTE) begin
                                    idx           <= 8'd0;
                                    reporter_busy <= 1'b0;
                                    state         <= ST_IDLE;
                                end
                                else begin
                                    idx   <= idx + 8'd1;
                                    state <= ST_SEND;
                                end
                            end
                        end

                        default: state <= ST_IDLE;
                    endcase
                end
            end
        end
    end

endmodule
