// ============================================================
// tb_cmd_fifo.v — verification of top_module's qubit command FIFO [v11]
// ============================================================
// THE PROBLEM UNDER TEST
//   v9/v10 held the qubit command in a SINGLE register. With the PC sending
//   back-to-back, bytes are 86.8 µs apart (10 bits @115200) while the FPGA
//   needs ≥220 µs per qubit ⇒ a byte arriving while the FSM is busy OVERWRITES
//   the previous one and disappears. Measured effect: P_click ~2.8× below the
//   model (bb84_uwoc_measure.py --chunk 64), sifted bits fall with it, QBER often 0%.
//
// THE TEST
//   Fire N_CMD qubit commands BACK-TO-BACK with no gap, then count how many
//   qubits the FPGA actually ran = chan_inst.click_count + chan_inst.lost_count
//   (both counters increment on every sample_pulse, i.e. each FSM_ENCODE entry).
//   PASS ⟺ qubits run == N_CMD  and  cmd_drop_cnt == 0.
//
//   vlog +incdir+. *.v && vsim -c -do "run -all; quit -f" tb_cmd_fifo
//
// RECORDED RESULTS (ModelSim ASE 10.1d, clear ocean d = 5 m, L5)
//   old RTL N_CMD = 32 → only 14/32 ran (56% lost, MATCHING the ~2.8× P_click gap)
//   new RTL N_CMD = 32 → 32/32,  drop 0
//           N_CMD = 64 → 64/64,  drop 0
//           N_CMD = 160→ 127/160, drop 33 — the FIFO overflows but it is COUNTED
//                        and LEDG[7] lights up; no more silent command loss.
//   (The FIFO drains while the PC is still sending, so the overflow threshold is
//    ~105 commands for a batch with no clicks; each click costs a further ~4 ms
//    ≈ 46 bytes, so 64 is the practical safe level.)
// ============================================================
`timescale 1ns / 1ps

// N_CMD corresponds to --chunk of bb84_uwoc_measure.py. Change it with (note the
// UPPERCASE -G; lowercase -g is mis-parsed by ModelSim 10.1d):
//     vsim -c -GN_CMD=64 -do "run -all; quit -f" tb_cmd_fifo
module tb_cmd_fifo #(
    parameter integer N_CMD = 32         // < 64 (FIFO depth) even when clicks occur
);

    // 50 MHz → 20 ns; 460800 baud → 108.5 cycles = 2170 ns / bit.
    // [Q2] Raised from 115200 (8680 ns) together with uart_rx/uart_tx in
    // top_module: at real time a 96-byte security packet costs 8.33 ms at
    // 115200 against a 6.55 ms block, so the UART would be the bottleneck.
    // This testbench models the HOST, so its bit period has to track the link.
    // Leaving it at 8680 makes uart_rx decode garbage, some of it with bit[7]
    // set, which lands in the qubit FIFO and the TB reports MORE qubits run
    // than commands sent.
    localparam integer BIT_NS = 2170;

    reg         clk = 1'b0;
    reg  [9:0]  sw;
    reg  [3:0]  key;
    reg         uart_rxd = 1'b1;
    wire [9:0]  ledr;
    wire [7:0]  ledg;
    wire [6:0]  hex0, hex1, hex2, hex3;
    wire        uart_txd;
    wire [35:0] gpio0, gpio1;

    always #10 clk = ~clk;

    top_module dut (
        .CLOCK_50(clk), .SW(sw), .KEY(key),
        .LEDR(ledr), .LEDG(ledg),
        .HEX0(hex0), .HEX1(hex1), .HEX2(hex2), .HEX3(hex3),
        .UART_TXD(uart_txd), .UART_RXD(uart_rxd),
        .GPIO_0(gpio0), .GPIO_1(gpio1)
    );

    // ---- Send one 8N1 byte, with NO idle gap between bytes ----
    task uart_send(input [7:0] b);
        integer i;
        begin
            uart_rxd = 1'b0;                     // start
            #(BIT_NS);
            for (i = 0; i < 8; i = i + 1) begin
                uart_rxd = b[i];                 // LSB first
                #(BIT_NS);
            end
            uart_rxd = 1'b1;                     // stop
            #(BIT_NS);
        end
    endtask

    integer k;
    integer executed;
    reg [7:0] cmd;

    initial begin
        // SW[9]=1 PC mode, SW[4]=1 channel on, SW[3:2]=00 clear ocean,
        // SW[7:5]=101 turbulence L5, SW[1]=0 adaptive off (μ, λ fixed)
        sw  = 10'b10_1011_0000;
        sw[4] = 1'b1;
        key = 4'b0111;                           // KEY[3]=rst_n=0
        #2000;
        key = 4'b1111;                           // release reset
        #20000;

        uart_send(8'h30);                        // dist_idx = 0  (d = 5 m)
        uart_send(8'h40);                        // water = 0, λ = 450 nm
        uart_send(8'h55);                        // turb = 5
        uart_send(8'h01);                        // reset statistics + flush the FIFO
        #200000;

        $display("[TB] ban %0d lenh qubit lien tiep, khong nghi...", N_CMD);
        for (k = 0; k < N_CMD; k = k + 1) begin
            // bit[7]=1 lenh qubit; bit[2]=a_data, bit[1]=a_basis, bit[0]=b_basis
            cmd = 8'h80 | (k[2:0]);
            uart_send(cmd);
        end

        // Wait for the FPGA to drain the FIFO: N_CMD × 220 µs, plus margin for the
        // clicks (each click costs a further ~4 ms in FSM_REPORT).
        #(N_CMD * 500_000 + 20_000_000);

        executed = dut.chan_inst.click_count + dut.chan_inst.lost_count;
        $display("[TB] lenh da gui        : %0d", N_CMD);
        $display("[TB] qubit da chay      : %0d  (click=%0d, lost=%0d)",
                 executed, dut.chan_inst.click_count, dut.chan_inst.lost_count);
        $display("[TB] cmd_drop_cnt       : %0d", dut.cmd_drop_cnt);
        $display("[TB] LEDG[7] (co ban)   : %0b", ledg[7]);

        if (executed == N_CMD && dut.cmd_drop_cnt == 0) begin
            $display("[TB] ==> DAT: khong mat lenh nao.");
        end else begin
            $display("[TB] ==> HONG: mat %0d lenh (thanh ghi don / FIFO tran).",
                     N_CMD - executed);
        end
        $finish;
    end

endmodule
