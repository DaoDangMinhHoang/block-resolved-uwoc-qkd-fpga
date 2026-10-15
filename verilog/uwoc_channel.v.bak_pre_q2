// ============================================================
// UWOC CHANNEL EMULATOR — Underwater Optical Wireless Channel
// Replaces gamma_gamma_channel (the atmospheric FSO channel)
// ============================================================
//
// PHYSICAL MODEL
//   h = h_l(d, λ; w) · h_f,  h_f = h_s · h_o      [Salcedo-Serrano ICC'22 eq.2]
//     h_l  : deterministic transmittance (Beer-Lambert + geometry) — psig_rom
//            w = water type; it and λ enter ONLY through c_w(λ)
//     h_s  : SCATTERING fading      ~ Gamma(1/σ_s², σ_s²) — hs_rom
//     h_o  : TURBULENCE fading      ~ Lognormal / Weibull — ho_rom
//     h_f  : fading pair, E[h_f] = 1  ⟹  E[h] = h_l  (h CONTAINS h_l, h_f does not)
//
// CORE DIFFERENCES from gamma_gamma_final.v
//   [1] There is a PHOTON DETECTION STAGE. Underwater most pulses give NO click
//       (P_click ~ 10⁻³…10⁻⁷), quite unlike the atmosphere where P_click ~ O(1).
//       → the dominant event is PHOTON LOSS, not "deep fade → random bit".
//   [2] There is a POLARISATION ERROR from multiple scattering e_pol(d) — the
//       UWOC-specific mechanism [1] Kebapci'23 blames for QBER rising with range.
//   [3] Fading is sampled on the COHERENCE TIME, NOT per qubit.
//       gamma_gamma_final.v resampled every qubit (sample_en) — physically wrong:
//       ocean turbulence varies ~10⁴× slower than the qubit period, so within
//       one monitoring window h must be treated as a constant. Resampling per
//       qubit averages the fading away and leaves adaptive_controller "blind".
//
//       ⚠ [v12] THE COHERENCE PERIOD IS NOW COUNTED IN QUBIT EVENTS, NOT CLOCK
//       CYCLES. What the physics fixes is the DIMENSIONLESS ratio
//           N_coh = τ_coh · f_rep = pulses per coherence block,
//       which the model puts at 5 ms × 10 MHz = 50 000 (LinkConfig.coherence_pulses).
//       A free-running clock timer only reproduces that when qubits really arrive
//       at f_rep. Driven from the PC over UART the board runs ~3500 qubit/s, so
//       2^18 cycles ≈ 5.24 ms held h constant for just 18 qubits — 2726× too few.
//       Every measurement window then spanned 10^5 independent fading samples and
//       averaged the turbulence away: measured QBER came out identical at L1 and
//       L5, and the variance the adaptive controller is supposed to track was not
//       observable at all.
//
//       coh_sel now selects the mode:
//           coh_sel = 0      legacy — resample every 2^COH_LOG2 CLOCK CYCLES
//           coh_sel = k ≥ 1  resample every 2^(k+5) QUBIT EVENTS
//       k = 11 gives 2^16 = 65 536 qubits, i.e. exactly one channel_monitor
//       window (NEXP_LOG2 = 16) per fading sample — the alignment that makes the
//       per-window QBER statistics comparable with the Python model called as
//       UWOCChannel.window_statistics(window = 2^(k+5)).
//   [4] h uses 12 bits (256 = 1.0), probabilities 24 bits — see uwoc_lut_gen.py
//       for why 8-bit/16-bit are NOT enough (heavy tails + tiny P_sig).
//
// DECISION STAGE (exactly equivalent to the QBER formula, NO division needed)
//       sig_det   = rand24 < p_sig        // click from the signal
//       noise_det = rand24 < p_noise      // click from dark count + background
//       click     = sig_det | noise_det
//       err       = sig_det ? (rand24 < e_pol) : rand24[23]   // noise → 50%
//   Yields exactly:  QBER = [e_pol·(1−e^(−n̄)) + ½Y₀] / P_click
//
// SIGNAL-PATH INTERFACE (OOK frame [SYNC=1][b1=basis][b0=data][IDLE=0])
//   • no_click  → force signal_out = 0 for the whole frame ⇒ the RX sees no
//                 rising edge ⇒ the top_module FSM times out ⇒ evt_qubit_lost.
//                 Exactly the semantics of "the photon never arrived".
//   • err       → flip the DATA BIT (slot b0) only. The SYNC slot is untouched
//                 (that would turn the frame into a lost photon) and so is the
//                 basis slot (a polarisation error mis-gates the detector WITHIN
//                 THE SAME basis — just like the spy model in bob.v).
//
// ESTIMATED RESOURCES: ~46 kbit ROM (11 M4K blocks) + 2 multipliers + ~400 LE
//
// REFERENCES
//   [1] Kebapci et al., IEEE Photonics J. 15(4), 2023 (FPGA underwater QKD)
//   [2] Salcedo-Serrano et al., IEEE ICC 2022, pp. 3814–3819
//   [3] Jamali, Akhoundi, Salehi, IEEE TWC 15(6), 2016
// ============================================================

module uwoc_channel #(
    // Fading resample period = the coherence time.
    // Default 2^18 cycles @50 MHz ≈ 5.24 ms — matches the model's τ_coh = 5 ms.
    parameter COH_LOG2 = 18
) (
    input  wire        clk,
    input  wire        rst_n,

    // ---- Signal path ----
    input  wire        signal_in,       // from the OOK TX (after PWM)
    output wire        signal_out,      // to the OOK RX
    input  wire [23:0] slot_width,      // used to locate the b0 data slot

    // ---- Channel configuration ----
    input  wire        chan_enable,     // 0 = ideal bypass (for testing)
    input  wire [1:0]  water_type,      // 0=clear ocean, 1=coastal, 2=harbor
    input  wire [3:0]  dist_idx,        // range index in the water type's grid
    input  wire [2:0]  turb_level,      // 0=OFF … 5=SEVERE
    input  wire [1:0]  lambda_sel,      // 0=450nm blue, 1=532nm green, 2=650nm red
    input  wire [3:0]  mu_level,        // transmit intensity; 8 = nominal μ
    input  wire [23:0] p_noise,         // Y₀ (dark+background), thang 2^24
    input  wire        dynamic_enable,  // allow the turbulence level to random-walk

    // ---- [v12] Coherence-block selector ----
    input  wire [3:0]  coh_sel,         // 0 = legacy 2^COH_LOG2 clock cycles;
                                        // k ≥ 1 = 2^(k+5) QUBIT EVENTS
    input  wire        stats_rst,       // pulse: the host reset its counters. Restart
                                        // the block grid here so the host can bin a
                                        // click by (total_qubits >> (coh_sel+5)).

    // ---- Status outputs ----
    output reg  [2:0]  current_level,   // effective turbulence level (range-compensated)
    output reg         click,           // pulse: this qubit DID click
    output reg         no_click,        // pulse: this qubit was LOST (no click)
    output reg         err_inject,      // this qubit has a data-bit error
    output reg  [15:0] click_count,
    output reg  [15:0] lost_count,
    output reg  [15:0] flip_count,

    // ---- Telemetry (for UART / channel_monitor) ----
    output reg  [11:0] h_s_out,         // 256 = 1.0
    output reg  [11:0] h_o_out,         // 256 = 1.0
    output reg  [11:0] h_f_out,         // h_s·h_o, 256 = 1.0
    output reg  [23:0] qub_index,       // [v12] 0-based ATTEMPT index of the qubit
                                        // just sampled, reset by stats_rst. The
                                        // top level's total_cnt cannot serve this
                                        // role: a lost photon skips FSM_PROCESS,
                                        // so total_cnt counts CLICKS, not attempts.
                                        // Coherence block = qub_index >> (coh_sel+5).
    output wire [15:0] nexp_inv,        // 32768 / E[clicks per window]

    // ---- Sample strobe, once per qubit (from the top-level FSM_ENCODE) ----
    input  wire        sample_en
);

    localparam H_MEAN = 12'd256;        // 256 ↔ h = 1.0
    localparam H_MAXV = 12'd4095;
    localparam [23:0] P_MAXV = 24'hFF_FFFF;

    // ============================================================
    // 1. RANDOM SOURCES
    // ============================================================
    // 5 independent LFSRs, long enough not to repeat within a measurement run.
    // (This is a channel emulator, NOT a cryptographic random source — the real
    //  entropy for BB84 still comes from trng.v.)

    reg [31:0] lfsr_hs, lfsr_ho;
    reg [31:0] lfsr_sig, lfsr_noise, lfsr_err;

    wire fb_hs    = lfsr_hs[31]    ^ lfsr_hs[21]    ^ lfsr_hs[1]    ^ lfsr_hs[0];
    wire fb_ho    = lfsr_ho[31]    ^ lfsr_ho[28]    ^ lfsr_ho[19]   ^ lfsr_ho[0];
    wire fb_sig   = lfsr_sig[31]   ^ lfsr_sig[30]   ^ lfsr_sig[10]  ^ lfsr_sig[0];
    wire fb_noise = lfsr_noise[31] ^ lfsr_noise[27] ^ lfsr_noise[5] ^ lfsr_noise[0];
    wire fb_err   = lfsr_err[31]   ^ lfsr_err[25]   ^ lfsr_err[13]  ^ lfsr_err[0];

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            lfsr_hs    <= 32'hDEAD_BEEF;
            lfsr_ho    <= 32'hCAFE_1337;
            lfsr_sig   <= 32'h1234_ABCD;
            lfsr_noise <= 32'h5A5A_F0F0;
            lfsr_err   <= 32'h0BAD_C0DE;
        end else begin
            lfsr_hs    <= {lfsr_hs[30:0],    fb_hs};
            lfsr_ho    <= {lfsr_ho[30:0],    fb_ho};
            lfsr_sig   <= {lfsr_sig[30:0],   fb_sig};
            lfsr_noise <= {lfsr_noise[30:0], fb_noise};
            lfsr_err   <= {lfsr_err[30:0],   fb_err};
        end
    end

    // ============================================================
    // 2. LOOKUP TABLES (generated by python/uwoc_lut_gen.py)
    // ============================================================
    `include "uwoc_channel_rom.vh"

    // The tables hold only 3 water types and 3 wavelengths; clamp so a value of 3 stays in range.
    wire [1:0] water_q  = (water_type > 2'd2) ? 2'd2 : water_type;
    wire [1:0] lam_q    = (lambda_sel > 2'd2) ? 2'd2 : lambda_sel;

    wire [5:0] wd_addr  = {water_q, dist_idx};              // 3×16   = 0..47
    // λ-dependent tables: flat index = (lam·3 + water)·16 + dist.
    // Multiply by 3 with add-and-shift to avoid inferring a multiplier.
    wire [3:0] lw_idx   = {2'b00, lam_q} + {1'b0, lam_q, 1'b0}; // lam·3
    wire [7:0] lwd_addr = {(lw_idx + {2'b00, water_q}), dist_idx};

    // ⚠ ROM READS MUST BE SYNCHRONOUS.
    //   Cyclone II M4K blocks only support registered reads. A combinational read
    //   (assigned straight to a wire) forces Quartus to unroll the tables into LEs:
    //   psig+epol+nexp_inv alone is 144 entries × 64 bits ≈ 9.2 kbit of combinational
    //   logic, an estimated 2000–3000 LE. Worse, once unrolled into logic the
    //   `initial` block content risks not being kept → ROM reads 0 → P_sig = 0 →
    //   the emulator NEVER clicks, and nothing reports an error.
    //   Registering these reads is functionally harmless: they only change when the
    //   channel configuration changes (rare, set by the user over UART/switches), so
    //   one extra cycle of latency is irrelevant. In return M4Ks are inferred and the
    //   ROM content follows Altera's standard memory-initialisation path.
    //   hs_rom / ho_rom need NO change — they are already read in a clocked block.
    reg [2:0]  hs_class;
    reg [23:0] psig_ref;
    reg [23:0] epol_p;
    reg [2:0]  ho_off_u;                           // 3-bit two's complement
    reg [15:0] nexp_inv_r;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            hs_class   <= 3'd0;
            psig_ref   <= 24'd0;
            epol_p     <= 24'd0;
            ho_off_u   <= 3'd0;
            nexp_inv_r <= 16'd32768;
        end else begin
            hs_class   <= hscls_rom[wd_addr];
            psig_ref   <= psig_rom[lwd_addr];
            epol_p     <= epol_rom[lwd_addr];
            ho_off_u   <= hooff_rom[wd_addr];
            nexp_inv_r <= nexp_inv_rom[lwd_addr];
        end
    end

    assign nexp_inv = nexp_inv_r;

    // Range-based turbulence offset: eff = clamp(turb_level + offset, 0, 5).
    // ho_rom is generated only at the 20 m reference range; the offset (log₄ of the
    // σ² ratio) maps the level onto σ²_ho(d) without a ROM per range.
    wire signed [4:0] ho_off  = {{2{ho_off_u[2]}}, ho_off_u};   // sign extension
    wire signed [5:0] lvl_raw = $signed({3'b000, turb_level}) + ho_off;

    // ============================================================
    // 3. EFFECTIVE TURBULENCE LEVEL (+ optional random walk)
    // ============================================================
    reg [25:0] dyn_timer;
    reg  [2:0] walk_level;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            dyn_timer  <= 26'd0;
            walk_level <= 3'd0;
        end else begin
            dyn_timer <= dyn_timer + 1'b1;
            if (!chan_enable) begin
                walk_level <= 3'd0;
            end else if (dyn_timer == 26'd0 && dynamic_enable) begin
                case (lfsr_ho[1:0])
                    2'b01:   walk_level <= (walk_level < 3'd1) ? walk_level + 1'b1 : 3'd1;
                    2'b10:   walk_level <= (walk_level > 3'd0) ? walk_level - 1'b1 : 3'd0;
                    default: walk_level <= walk_level;
                endcase
            end else if (!dynamic_enable) begin
                walk_level <= 3'd0;
            end
        end
    end

    // ⚠ The extension MUST be SIGNED. Using {1'b0, lvl_raw} zero-extends a
    //   negative value (e.g. −1 → +31) and pins the turbulence level at the ceiling
    //   — a bug the testbench caught (eff_lv = 5 where it should have been 1).
    wire signed [6:0] walk_ext = $signed({4'b0000, walk_level});
    wire signed [6:0] lvl_sum  = lvl_raw + walk_ext;

    always @(*) begin
        if (!chan_enable || turb_level == 3'd0)
            current_level = 3'd0;
        else if (lvl_sum < 7'sd1)
            current_level = 3'd1;
        else if (lvl_sum > 7'sd5)
            current_level = 3'd5;
        else
            current_level = lvl_sum[2:0];
    end

    // ============================================================
    // 4. FADING SAMPLED ON THE COHERENCE TIME
    // ============================================================
    // ⚠ THIS IS THE MOST IMPORTANT DIFFERENCE from gamma_gamma_final.v.
    //   The old module resampled h every qubit (sample_en) → fading was averaged
    //   out over every monitoring window → the controller saw no turbulence.
    //   Here h_s and h_o stay frozen for a whole coherence block.
    //
    //   [v12] The block length is counted in QUBIT EVENTS when coh_sel ≠ 0 — see
    //   the note [3] in the header for why the clock-cycle version silently
    //   destroyed the turbulence signal whenever the qubit rate was below f_rep.
    //   coh_sel = 0 keeps the old clock timer so previously collected data sets
    //   remain reproducible bit for bit.

    // Declared here rather than in section 5: the coherence counter below needs
    // the once-per-qubit strobe, and Verilog-2001 wants the net in scope first.
    reg  sample_en_d;
    wire sample_pulse = sample_en & ~sample_en_d;   // rising edge only

    reg [COH_LOG2-1:0] coh_timer;
    reg [23:0]         qub_cnt;      // qubit events since the last stats reset
    reg                coh_bit_d;
    reg [11:0] h_s_reg, h_o_reg;

    // Bit n of a free-running counter toggles every 2^n increments, so watching
    // qub_cnt[coh_sel+5] gives a period of 2^(coh_sel+5) qubits for the price of a
    // 24:1 mux — no wide comparator, no divider.
    //   coh_sel = 1 →     64 qubits      coh_sel = 11 →  65 536 qubits
    //   coh_sel = 8 →  8 192 qubits      coh_sel = 15 → 1 048 576 qubits
    wire [4:0] coh_log2_q = {1'b0, coh_sel} + 5'd5;
    wire       coh_bit    = qub_cnt[coh_log2_q];

    wire coh_tick_clk = (coh_sel == 4'd0) && (coh_timer == {COH_LOG2{1'b0}});
    wire coh_tick_qub = (coh_sel != 4'd0) && (coh_bit != coh_bit_d);
    wire coh_tick     = coh_tick_clk | coh_tick_qub;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            coh_timer <= {COH_LOG2{1'b0}};
            qub_cnt   <= 24'd0;
            coh_bit_d <= 1'b0;
            h_s_reg   <= H_MEAN;
            h_o_reg   <= H_MEAN;
        end else if (!chan_enable) begin
            coh_timer <= {COH_LOG2{1'b0}};
            qub_cnt   <= 24'd0;
            coh_bit_d <= 1'b0;
            h_s_reg   <= H_MEAN;
            h_o_reg   <= H_MEAN;
        end else if (stats_rst) begin
            // Draw a FRESH h here too. Without it the first block of a new
            // measurement point would inherit the fading sample of the previous
            // configuration — which at 2^16 qubits per block is no longer a
            // negligible fraction of the point.
            coh_timer <= {COH_LOG2{1'b0}};
            qub_cnt   <= 24'd0;
            coh_bit_d <= 1'b0;
            h_s_reg   <= hs_rom[{hs_class,      lfsr_hs[7:0]}];
            h_o_reg   <= ho_rom[{current_level, lfsr_ho[7:0]}];
        end else begin
            coh_timer <= coh_timer + 1'b1;
            if (sample_pulse)
                qub_cnt <= qub_cnt + 24'd1;
            coh_bit_d <= coh_bit;
            if (coh_tick) begin
                h_s_reg <= hs_rom[{hs_class,      lfsr_hs[7:0]}];
                h_o_reg <= ho_rom[{current_level, lfsr_ho[7:0]}];
            end
        end
    end

    // ---- Step 1: h_f = (h_s · h_o) >> 8, saturated to 12 bits ----
    wire [23:0] hf_prod = h_s_reg * h_o_reg;            // 12×12 → 24
    wire [15:0] hf_shift = hf_prod[23:8];
    wire [11:0] h_f = (hf_shift > {4'd0, H_MAXV}) ? H_MAXV : hf_shift[11:0];

    // ---- Step 2: p_sig = ((psig_ref · h_f) >> 8) · mu_level >> 3 ----
    // Linearising P_sig ≈ μ·h_l·h_f·η is valid because n̄ ≪ 1 underwater
    // (check 6a in uwoc_channel_model.py: the Jensen gap is < 1%).
    wire [35:0] ps_prod  = psig_ref * h_f;              // 24×12 → 36
    wire [27:0] ps_shift = ps_prod[35:8];
    wire [23:0] ps_h     = (ps_shift[27:24] != 4'd0) ? P_MAXV : ps_shift[23:0];

    wire [27:0] ps_mu    = ps_h * mu_level;             // 24×4 → 28
    wire [24:0] ps_mu_s  = ps_mu[27:3];
    wire [23:0] p_sig    = (ps_mu_s[24]) ? P_MAXV : ps_mu_s[23:0];

    // ============================================================
    // 5. CLICK / ERROR DECISION — ONCE PER QUBIT
    // ============================================================
    // sample_en_d / sample_pulse are declared in section 4 (the coherence counter
    // needs them); sample_en_d is still registered in the always block below.
    wire [23:0] rnd_sig   = lfsr_sig[23:0];
    wire [23:0] rnd_noise = lfsr_noise[23:0];
    wire [23:0] rnd_err   = lfsr_err[23:0];

    wire sig_det_w   = (rnd_sig   < p_sig);
    wire noise_det_w = (rnd_noise < p_noise);
    wire click_w     = sig_det_w | noise_det_w;
    // A signal-induced click → wrong with probability e_pol.
    // A background-only click → fully random ⇒ wrong 50% of the time (1 bit).
    wire err_w       = sig_det_w ? (rnd_err < epol_p) : rnd_err[23];

    reg qubit_click, qubit_err;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sample_en_d <= 1'b0;
            qubit_click <= 1'b1;
            qubit_err   <= 1'b0;
            click       <= 1'b0;
            no_click    <= 1'b0;
            err_inject  <= 1'b0;
            click_count <= 16'd0;
            lost_count  <= 16'd0;
            flip_count  <= 16'd0;
            h_s_out     <= H_MEAN;
            h_o_out     <= H_MEAN;
            h_f_out     <= H_MEAN;
            qub_index   <= 24'd0;
        end else begin
            sample_en_d <= sample_en;
            click       <= 1'b0;
            no_click    <= 1'b0;
            err_inject  <= 1'b0;

            if (sample_pulse) begin
                if (!chan_enable) begin
                    // Bypass: ideal channel (used to test everything else)
                    qubit_click <= 1'b1;
                    qubit_err   <= 1'b0;
                    click       <= 1'b1;
                end else begin
                    qubit_click <= click_w;
                    qubit_err   <= err_w & click_w;
                    click       <=  click_w;
                    no_click    <= ~click_w;
                    err_inject  <=  err_w & click_w;

                    if (click_w) begin
                        click_count <= click_count + 1'b1;
                        if (err_w) flip_count <= flip_count + 1'b1;
                    end else begin
                        lost_count <= lost_count + 1'b1;
                    end
                end
                h_s_out <= h_s_reg;
                h_o_out <= h_o_reg;
                h_f_out <= h_f;
                // PRE-increment value: qub_cnt is bumped by the coherence block
                // below at this same edge, so this latches the 0-based index of
                // the qubit being sampled right now. It therefore satisfies
                //     block(this qubit) == qub_index >> (coh_sel + 5)
                // exactly, with no off-by-one against the h that was applied.
                qub_index <= qub_cnt;
            end
        end
    end

    // ============================================================
    // 6. APPLY TO THE SIGNAL PATH
    // ============================================================
    // Tracks the OOK frame: [SYNC][b1=basis][b0=data][IDLE], each slot = slot_width.
    // Counting starts on the SYNC rising edge; the flip window is slot 3 (b0), i.e.
    // t ∈ [2·slot_width, 3·slot_width).

    reg  sig_prev;
    reg  [25:0] frame_cnt;
    reg  frame_run;

    wire [25:0] t_b0_start = {slot_width[23:0], 1'b0};              // 2·W
    wire [25:0] t_b0_end   = {slot_width[23:0], 1'b0} + {2'b00, slot_width}; // 3·W

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sig_prev  <= 1'b0;
            frame_cnt <= 26'd0;
            frame_run <= 1'b0;
        end else begin
            sig_prev <= signal_in;
            if (signal_in & ~sig_prev & ~frame_run) begin
                frame_cnt <= 26'd0;
                frame_run <= 1'b1;
            end else if (frame_run) begin
                if (frame_cnt >= t_b0_end)
                    frame_run <= 1'b0;
                else
                    frame_cnt <= frame_cnt + 1'b1;
            end
        end
    end

    wire in_data_slot = frame_run
                      && (frame_cnt >= t_b0_start)
                      && (frame_cnt <  t_b0_end);

    // no_click → silence for the whole frame (the RX times out ⇒ evt_qubit_lost)
    // click + error → flip the data slot only
    assign signal_out = (!chan_enable)  ? signal_in
                      : (!qubit_click)  ? 1'b0
                      : (qubit_err && in_data_slot) ? ~signal_in
                                                    :  signal_in;

endmodule
