#!/bin/bash
# For the reaction field case
# Run Drude-Lorentz reference calculation:
cd RF_drl_ref
../../../WaveT/bin/WaveT-serial.x < tdcis.inp > out.dat
# Run Drude-Lorentz calculation with the algorithm for general dielectric functions:
# the file eps.inp contains the input complex dielectric function table
cd ../RF_drl_gen
../../../WaveT/bin/WaveT-serial.x < tdcis.inp > out.dat
# For the local field case
# Run Drude-Lorentz reference calculation:
cd ../LF_drl_ref
../../../WaveT/bin/WaveT-serial.x < tdcis.inp > out.dat
# Run Drude-Lorentz calculation with the algorithm for general dielectric functions:
# the file eps.inp contains the input complex dielectric function table
cd ../LF_drl_gen
../../../WaveT/bin/WaveT-serial.x < tdcis.inp > out.dat
cd ../
#Compare results:
gnuplot plot.gnu

