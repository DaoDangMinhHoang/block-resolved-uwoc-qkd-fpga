# ==============================================================================
# ModelSim DO script - UWOC Channel Emulator Testbenches
# Dung trong ModelSim Console: do run_sim.do [tb_name]
# Vi du:
#   do run_sim.do tb_security_chain
#   do run_sim.do tb_index_align
#   do run_sim.do tb_uwoc_channel
#   do run_sim.do tb_cmd_fifo
#   do run_sim.do tb_adaptive_loop
#   do run_sim.do tb_realtime_attempt_engine
# ==============================================================================

if {![file exists work]} {
    vlib work
}

vlog +incdir+. *.v

if {$argc > 0} {
    set tb [lindex $argv 0]
} else {
    set tb "tb_security_chain"
}

puts ">> Dang khoi chay mo phong: $tb"
vsim work.$tb
add wave -r /*
run -all
