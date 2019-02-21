#!/bin/bash
# Run equation of motion stemming directly from Onsager model (of a point-dipole inside an spherical cavity)
cd DIP-Sphere
../../../bin/WaveT-serial.x < tdcis.inp > out.dat
# Run CPCM-like polarization charges equation motion with Onsager relaxation time and considering the interaction with just a molecular dipole
cd ../ONS-ONS
#cp out_write/* .
../../../bin/WaveT-serial.x < tdcis.inp > out.dat
# Run CPCM-like polarization charges equation motion with Onsager relaxation time and considering the interaction with complete molecular charge distribution
cd ../ONS-PCM
#cp out_write/* .
../../../bin/WaveT-serial.x < tdcis.inp > out.dat
# Run IEF-PCM polarization charges equation motion considering the interaction with a molecular dipole
cd ../IEF-ONS
#cp out_write/* .
../../../bin/WaveT-serial.x < tdcis.inp > out.dat
# Run IEF-PCM polarization charges equation motion considering the interaction with complete molecular charge distribution
cd ../IEF-PCM
../../../bin/WaveT-serial.x < tdcis.inp > out.dat
cd ../
#Visualize results:
#1) Test results:
gnuplot plot.gnu
#2) Compare test results with references:
gnuplot compare.gnu

