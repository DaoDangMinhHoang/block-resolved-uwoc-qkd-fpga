# Pin assignments for Terasic DE1 Board (Cyclone II EP2C20F484C7)
# For project: QKD FSO (top_module)

# Clock
set_location_assignment PIN_L1 -to CLOCK_50

# Switches (SW)
set_location_assignment PIN_L22 -to SW[0]
set_location_assignment PIN_L21 -to SW[1]
set_location_assignment PIN_M22 -to SW[2]
set_location_assignment PIN_V12 -to SW[3]
set_location_assignment PIN_W12 -to SW[4]
set_location_assignment PIN_U12 -to SW[5]
set_location_assignment PIN_U11 -to SW[6]
set_location_assignment PIN_M2 -to SW[7]
set_location_assignment PIN_M1 -to SW[8]
set_location_assignment PIN_L2 -to SW[9]

# Keys / Pushbuttons (Active Low)
set_location_assignment PIN_R22 -to KEY[0]
set_location_assignment PIN_R21 -to KEY[1]
set_location_assignment PIN_T22 -to KEY[2]
set_location_assignment PIN_T21 -to KEY[3]

# Red LEDs
set_location_assignment PIN_R20 -to LEDR[0]
set_location_assignment PIN_R19 -to LEDR[1]
set_location_assignment PIN_U19 -to LEDR[2]
set_location_assignment PIN_Y19 -to LEDR[3]
set_location_assignment PIN_T18 -to LEDR[4]
set_location_assignment PIN_V19 -to LEDR[5]
set_location_assignment PIN_Y18 -to LEDR[6]
set_location_assignment PIN_U18 -to LEDR[7]
set_location_assignment PIN_R18 -to LEDR[8]
set_location_assignment PIN_R17 -to LEDR[9]

# Green LEDs
set_location_assignment PIN_U22 -to LEDG[0]
set_location_assignment PIN_U21 -to LEDG[1]
set_location_assignment PIN_V22 -to LEDG[2]
set_location_assignment PIN_V21 -to LEDG[3]
set_location_assignment PIN_W22 -to LEDG[4]
set_location_assignment PIN_W21 -to LEDG[5]
set_location_assignment PIN_Y22 -to LEDG[6]
set_location_assignment PIN_Y21 -to LEDG[7]

# 7-Segment Display 0
set_location_assignment PIN_J2 -to HEX0[0]
set_location_assignment PIN_J1 -to HEX0[1]
set_location_assignment PIN_H2 -to HEX0[2]
set_location_assignment PIN_H1 -to HEX0[3]
set_location_assignment PIN_F2 -to HEX0[4]
set_location_assignment PIN_F1 -to HEX0[5]
set_location_assignment PIN_E2 -to HEX0[6]

# 7-Segment Display 1
set_location_assignment PIN_E1 -to HEX1[0]
set_location_assignment PIN_H6 -to HEX1[1]
set_location_assignment PIN_H5 -to HEX1[2]
set_location_assignment PIN_H4 -to HEX1[3]
set_location_assignment PIN_G3 -to HEX1[4]
set_location_assignment PIN_D2 -to HEX1[5]
set_location_assignment PIN_D1 -to HEX1[6]

# 7-Segment Display 2
set_location_assignment PIN_G5 -to HEX2[0]
set_location_assignment PIN_G6 -to HEX2[1]
set_location_assignment PIN_C2 -to HEX2[2]
set_location_assignment PIN_C1 -to HEX2[3]
set_location_assignment PIN_E3 -to HEX2[4]
set_location_assignment PIN_E4 -to HEX2[5]
set_location_assignment PIN_D3 -to HEX2[6]

# 7-Segment Display 3
set_location_assignment PIN_F4 -to HEX3[0]
set_location_assignment PIN_D5 -to HEX3[1]
set_location_assignment PIN_D6 -to HEX3[2]
set_location_assignment PIN_J4 -to HEX3[3]
set_location_assignment PIN_L8 -to HEX3[4]
set_location_assignment PIN_F3 -to HEX3[5]
set_location_assignment PIN_D4 -to HEX3[6]

# UART (RS-232)
set_location_assignment PIN_F14 -to UART_RXD
set_location_assignment PIN_G12 -to UART_TXD

# GPIO 0 (used for loopback / optical signals)
set_location_assignment PIN_A13 -to GPIO_0[0]
set_location_assignment PIN_B13 -to GPIO_0[1]
