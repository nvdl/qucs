<Qucs Schematic 0.0.17>
<Properties>
  <View=90,-9,1329,668,0.832123,0,0>
  <Grid=10,10,1>
  <DataSet=X5_500.dat>
  <DataDisplay=X5_500.dpl>
  <OpenDisplay=1>
  <Script=X5_500.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <GND * 1 380 130 0 0 0 0>
  <GND * 1 610 80 0 0 0 0>
  <GND * 1 220 330 0 0 0 0>
  <GND * 1 510 330 0 0 0 0>
  <.TR TR1 1 650 50 0 68 0 0 "lin" 0 "0" 0 "1e-6" 1 "1001" 0 "" 0 "2" 0 "1 ns" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "0" 0>
  <GND * 1 380 400 0 0 0 0>
  <GND * 1 300 400 0 0 0 0>
  <R R5 5 300 370 -63 -9 0 1 "360" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <C C3 5 350 240 -46 -35 0 0 ".01 uF" 1 "" 0 "neutral" 0>
  <R R7 5 510 300 -84 -8 0 1 "3 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <R R1 5 440 170 -20 13 0 0 "2k" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <R R2 5 440 60 -20 -29 0 0 "100" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <L L4 5 280 240 -31 6 0 0 "58 nH" 1 "" 0>
  <_BJT BFR92A_1 5 510 240 -31 6 0 0 "pnp" 0 "4.11877E-016" 0 "9.97275E-001" 0 "9.96202E-001" 0 "3.20054E+000" 0 "1.28155E+000" 0 "6.26719E+001" 0 "3.36915E+000" 0 " 4.01062E-015" 0 "1.57708E+000" 0 "2.79905E-016" 0 "1.07543E+000" 0 "314" 0 "1.81086E+001" 0 "1.00000E+001" 0 "1.00000E-006" 0 "2.32000E+000" 0 "1.16450E+000" 0 "1.00000E+001" 0 "8.90512E-013" 0 "6.00000E-001" 0 "2.58570E-001" 0 "5.46563E-013" 0 "3.80824E-001" 0 "2.02935E-001" 0 "1" 0 "0" 0 "0.75" 0 "0" 0 "0.5" 0 "1.54973E-011" 0 "3.91402E+001" 0 "2.15279E+000" 0 "2.13776E-001" 0 "3.04e-07" 0 "26.85" 0 "0" 0 "1" 0 "1" 0 "0" 0 "1" 0 "1" 0 "0" 0 "1.5" 0 "3.00000E+000" 0 "1.11" 0 "26.85" 0 "1" 0>
  <Vdc V1 5 380 100 -71 -26 0 3 "15 V" 1>
  <.DC DC1 1 200 150 0 41 0 0 "26.85" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "no" 0 "150" 0 "no" 0 "none" 0 "CroutLU" 0>
  <C C4 5 570 170 -26 17 0 0 "0.5 pF" 1 "" 0 "neutral" 0>
  <GND * 1 870 260 0 0 0 0>
  <Pac P2 5 870 230 21 -21 0 1 "2" 0 "50 Ohm" 0 "-20 dBm" 0 "1 GHz" 0 "26.85" 0>
  <Eqn Eqn1 1 480 400 -30 16 0 0 "Spectrum_output_dBm=15+dB(Time2Freq(Vrf.Vt,time))" 1 "Spectrum_input_dBm=15+dB(Time2Freq(Vinput.Vt,time))" 1 "yes" 0>
  <R R8 5 810 170 -15 10 0 2 "0.3" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <C C1 5 560 60 -26 17 0 0 "56 pF" 1 "" 0 "neutral" 0>
  <L L1 5 510 110 -63 -6 0 3 "88 nH" 1 "" 0>
  <L L3 5 380 290 -75 -11 0 1 "24 nH" 1 "" 0>
  <C C2 5 380 370 21 -8 0 1 "3.3 pF" 1 "" 0 "neutral" 0>
  <L L6 5 730 170 -39 9 0 0 "228 nH" 1 "" 0>
  <L L5 5 610 260 -88 -12 0 1 "35.6 nH" 1 "" 0>
  <R R6 5 610 320 -60 -8 0 1 "0.1" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <C C6 5 690 300 10 22 0 1 "2.7 pF" 1 "" 0 "neutral" 0>
  <GND * 1 650 390 0 0 0 0>
  <Pac P1 5 220 290 -88 -7 0 1 "1" 0 "50 Ohm" 0 "0 dBm" 1 "100 MHz" 1 "27" 0>
</Components>
<Wires>
  <380 240 480 240 "" 0 0 0 "">
  <510 170 510 210 "" 0 0 0 "">
  <510 60 530 60 "" 0 0 0 "">
  <470 60 510 60 "" 0 0 0 "">
  <380 60 410 60 "" 0 0 0 "">
  <380 60 380 70 "" 0 0 0 "">
  <590 60 610 60 "" 0 0 0 "">
  <610 60 610 80 "" 0 0 0 "">
  <220 320 220 330 "" 0 0 0 "">
  <470 170 510 170 "" 0 0 0 "">
  <380 170 380 240 "" 0 0 0 "">
  <380 170 410 170 "" 0 0 0 "">
  <380 330 380 340 "" 0 0 0 "">
  <300 330 380 330 "" 0 0 0 "">
  <300 330 300 340 "" 0 0 0 "">
  <310 240 320 240 "" 0 0 0 "">
  <220 240 220 260 "" 0 0 0 "">
  <220 240 250 240 "" 0 0 0 "">
  <510 170 540 170 "" 0 0 0 "">
  <650 170 650 220 "" 0 0 0 "">
  <600 170 650 170 "" 0 0 0 "">
  <650 170 700 170 "" 0 0 0 "">
  <760 170 780 170 "" 0 0 0 "">
  <870 170 870 200 "" 0 0 0 "">
  <840 170 870 170 "" 0 0 0 "">
  <510 140 510 170 "" 0 0 0 "">
  <510 60 510 80 "" 0 0 0 "">
  <380 240 380 260 "" 0 0 0 "">
  <380 320 380 330 "" 0 0 0 "">
  <610 220 650 220 "" 0 0 0 "">
  <610 220 610 230 "" 0 0 0 "">
  <650 220 690 220 "" 0 0 0 "">
  <690 220 690 270 "" 0 0 0 "">
  <610 350 610 370 "" 0 0 0 "">
  <610 370 650 370 "" 0 0 0 "">
  <690 330 690 370 "" 0 0 0 "">
  <650 370 690 370 "" 0 0 0 "">
  <650 370 650 390 "" 0 0 0 "">
  <220 240 220 240 "Vinput" 130 170 0 "">
  <870 170 870 170 "Vrf" 780 250 0 "">
</Wires>
<Diagrams>
  <Rect 950 417 332 337 3 #c0c0c0 1 00 1 0 2e+08 1e+09 0 -50 10 10 1 -1 1 1 315 0 225 "Frequency, Hz" "" "">
	<"Spectrum_output_dBm" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
  <Text 530 330 12 #000000 0 "Inductor Q">
  <Text 780 140 12 #000000 0 "Inductor Q">
  <Text 190 490 14 #ff0000 0 "X 5 frequency multiplier implemented with a bipolar transistor, followed by a single pole filter. Although a lumped-constant\nfilter is shown, a high-Q coaxial or helical resonator structure would be more appropriate (followed perhaps by a 3 dB atten-\nuator and a conventional multi-pole bandpass filter to reduce harmonics to desired levels).  Note that series resistors approximate\ninductor "Q" for circuit loss/rejection estimates.  The circuit is adjusted by slightly bending coils of hand-wound air core inductors.">
  <Text 900 40 12 #000000 0 "Note ~ 5 dB gain from 100 MHz input to 500 MHz output">
</Paintings>