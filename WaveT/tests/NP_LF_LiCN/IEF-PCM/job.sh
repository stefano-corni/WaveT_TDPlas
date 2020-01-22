#!/bin/sh
#PBS -N pes_s1      
#PBS -j oe
#PBS -q s3par
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -W  group_list=moe-sc
#PBS -l mem=24000MB
cd  $PBS_O_WORKDIR

WaveT-serial.x < tdcis.inp 
