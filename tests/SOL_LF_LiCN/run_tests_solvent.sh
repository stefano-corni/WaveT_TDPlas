#!/bin/bash
# Run Onsager calculation:
cd DIP                 
../../../WaveT/bin/WaveT.x < tdcis.inp > out.dat
# Run Onsager calculation:
cd ../ONS
../../../WaveT/bin/WaveT.x < tdcis.inp > out.dat
# Run PCM calculation with only dipolar contributions (potentials calculated from dipoles)
cd ../PCM-mu
#cp out_write/* .
../../../WaveT/bin/WaveT.x < tdcis.inp > out.dat
# Run full PCM calculation
cd ../PCM
#cp out_write/* .
../../../WaveT/bin/WaveT.x < tdcis.inp > out.dat
cd ../
#Visualize results:
#1) Test results:
#gnuplot plot.gnu
#2) Compare test results with references:
#gnuplot compare.gnu
    

