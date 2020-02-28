#!/bin/sh
#PBS -N pes_s1      
#PBS -j oe
#PBS -q s3par
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -W  group_list=moe-sc
#PBS -l mem=24000MB
cd  $PBS_O_WORKDIR
./create_mesh.sh
tdplas.x < tdplas.inp
mv cavity.inp mat_SD.inp ../../NP_print_charges/print_charges/
cd ../../NP_print_charges/print_charges/
WaveT_serial.x < wavet.inp 
