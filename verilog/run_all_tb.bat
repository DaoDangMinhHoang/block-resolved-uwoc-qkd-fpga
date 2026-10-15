@echo off
setlocal enabledelayedexpansion

echo ==============================================================================
echo   UWOC FPGA Emulator - Chay Testbench tren ModelSim
echo ==============================================================================

cd /d "%~dp0"

if not exist work (
    echo [*] Tao thu vien work...
    vlib work
)

echo [*] Bien dich tat ca ma Verilog...
vlog +incdir+. *.v
if errorlevel 1 (
    echo [ERROR] Bien dich that bai!
    pause
    exit /b 1
)

if "%1"=="" goto run_all

echo [*] Chay testbench: %1
vsim -c -do "run -all; quit -f" %1
goto end

:run_all
echo.
echo ==============================================================================
echo [1/6] Chay tb_security_chain...
vsim -c -do "run -all; quit -f" tb_security_chain

echo.
echo ==============================================================================
echo [2/6] Chay tb_index_align...
vsim -c -do "run -all; quit -f" tb_index_align

echo.
echo ==============================================================================
echo [3/6] Chay tb_realtime_attempt_engine...
vsim -c -do "run -all; quit -f" tb_realtime_attempt_engine

echo.
echo ==============================================================================
echo [4/6] Chay tb_cmd_fifo...
vsim -c -do "run -all; quit -f" tb_cmd_fifo

echo.
echo ==============================================================================
echo [5/6] Chay tb_adaptive_loop...
vsim -c -do "run -all; quit -f" tb_adaptive_loop

echo.
echo ==============================================================================
echo [6/6] Chay tb_uwoc_channel...
vsim -c -do "run -all; quit -f" tb_uwoc_channel

:end
echo.
echo ==============================================================================
echo   HOAN TAT TAT CA TESTBENCH!
echo ==============================================================================
