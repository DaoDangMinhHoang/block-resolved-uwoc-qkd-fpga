# ============================================================
# bb84_phase2.sdc — Timing Constraints for BB84 QKD System
# Target: Altera Cyclone II EP2C20F484C7 (Terasic DE1)
# ============================================================
#
# FIX: the ring-oscillator TRNG creates deliberate combinational loops.
# TimeQuest cannot tell a deliberate loop from a design error
# → it reports negative setup slack. This file declares false paths for the ROs
# → TimeQuest then ignores those paths when computing timing.
#
# HOW TO ADD IT TO THE PROJECT:
# 1. Copy this file into the project directory (alongside the .qpf)
# 2. Quartus → Assignments → Settings → TimeQuest Timing Analyzer
#    → SDC files → Add → select bb84_phase2.sdc
# 3. Recompile (Processing → Start Compilation)
# 4. Open TimeQuest → verify setup slack > 0
# ============================================================

# ---- 1. Clock Definition ----
# DE1 board 50 MHz oscillator on PIN_L1
create_clock -name CLOCK_50 -period 20.000 [get_ports {CLOCK_50}]

# ---- 2. Clock Uncertainty ----
# Account for PLL jitter and board-level clock skew
derive_clock_uncertainty

# ---- 3. FALSE PATHS — Ring Oscillator TRNG ----
# Ring oscillators are intentional combinational loops used
# for entropy generation. They must be excluded from STA.
#
# 3 instances of trng_random → each contains trng with 4 ROs:
#   gen_alice_data/trng_core/ro{0,1,2,3}_chain[*]
#   gen_alice_basis/trng_core/ro{0,1,2,3}_chain[*]
#   gen_bob_basis/trng_core/ro{0,1,2,3}_chain[*]

# Method 1: Cut the asynchronous capture of the mixed RO output.
# This is the path that actually matters: raw_xor is a free-running
# combinational signal sampled by CLOCK_50 in trng.v.
#
# NOTE: do NOT use get_keepers for the RO chains themselves. get_keepers
# matches registers and ports only, while ro*_chain[*] are combinational
# LCELL nodes -> the collection comes back empty and TimeQuest silently
# drops the constraint (Warning 332174 / 332049).
set_false_path -to [get_registers {*trng_core*raw_sampled*}]

# Method 2: Cut every path leaving the ring oscillator chains.
# get_nodes matches combinational nodes, which is what the ROs are.
set_false_path -from [get_nodes {*trng_core*ro*_chain*}]

# ---- 4. FALSE PATHS — Asynchronous Inputs ----
# DIP switches and push buttons are asynchronous
set_false_path -from [get_ports {SW[*]}]
set_false_path -from [get_ports {KEY[*]}]

# UART RX input is asynchronous (oversampled by uart_rx module)
set_false_path -from [get_ports {UART_RXD}]

# ---- 5. Output Constraints ----
# LEDs and 7-segment have no critical timing
set_false_path -to [get_ports {LEDR[*]}]
set_false_path -to [get_ports {LEDG[*]}]
set_false_path -to [get_ports {HEX0[*]}]
set_false_path -to [get_ports {HEX1[*]}]
set_false_path -to [get_ports {HEX2[*]}]
set_false_path -to [get_ports {HEX3[*]}]

# UART TX output
set_false_path -to [get_ports {UART_TXD}]

# GPIO (active only in external laser mode)
set_false_path -to [get_ports {GPIO_0[*]}]
set_false_path -from [get_ports {GPIO_0[*]}]
set_false_path -to [get_ports {GPIO_1[*]}]
set_false_path -from [get_ports {GPIO_1[*]}]
