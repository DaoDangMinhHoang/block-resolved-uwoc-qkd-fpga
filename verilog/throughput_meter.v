// ============================================================================
// throughput_meter.v
// Counts attempt events over exactly CLK_HZ clocks and latches the total.
// Use this, not a Python wall-clock timer: UART and OS buffering make a
// host-side estimate meaningless at these rates.
//
// This is the measurement that proves the real-time claim. With
// realtime_attempt_engine paced at rate_hz = 10_000_000 the reported value must
// read 10_000_000 +/- 1. If it reads low, the channel cannot sustain the pace
// and engine_overruns will be nonzero.
//
// CLK_HZ MUST equal the real clock driving clk. If timing forces a lower clock,
// change it here and in the engine, or every reported rate is wrong by the
// same ratio.
// ============================================================================

module throughput_meter #(
    parameter integer CLK_HZ = 50000000
)(
    input  wire        clk,
    input  wire        rst_n,
    input  wire        clear,
    input  wire        evt_attempt,

    output reg  [31:0] attempts_per_sec,
    output reg         rate_valid,
    output reg         rate_pulse
);

    reg [31:0] clk_count;
    reg [31:0] attempt_count;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            clk_count        <= 32'd0;
            attempt_count    <= 32'd0;
            attempts_per_sec <= 32'd0;
            rate_valid       <= 1'b0;
            rate_pulse       <= 1'b0;
        end
        else if (clear) begin
            clk_count     <= 32'd0;
            attempt_count <= 32'd0;
            rate_pulse    <= 1'b0;
            // Do NOT retain the previous value. Rate presets now change between
            // campaigns, so a stale figure is wrong metadata: for the first
            // second, about 153 blocks at 10 Mattempt/s, every packet would
            // carry the rate of the run before. rate_valid tells the host to
            // ignore the field until the first full window has closed.
            attempts_per_sec <= 32'd0;
            rate_valid       <= 1'b0;
        end
        else begin
            rate_pulse <= 1'b0;

            if (evt_attempt)
                attempt_count <= attempt_count + 32'd1;

            if (clk_count == CLK_HZ - 1) begin
                attempts_per_sec <= attempt_count +
                                    (evt_attempt ? 32'd1 : 32'd0);
                rate_valid       <= 1'b1;
                attempt_count    <= 32'd0;
                clk_count        <= 32'd0;
                rate_pulse       <= 1'b1;
            end
            else begin
                clk_count <= clk_count + 32'd1;
            end
        end
    end

endmodule
