#!/bin/sh
#PBS -V
#PBS -N MPI 
#PBS -l nodes=1:ppn=4
#PBS -j oe
# starting stuff ...
set -x
date
env

export OMP_NUM_THREADS=1

# define a proper scratch temp dir
#
# put here commands
#
HOMEDIR=/home/coccia/
STORAGE=/storage/coccia
ROOT=OMP_tests/LiCN_NP/1_OMP



#TMPDIR=$STORAGE/scratch/$ROOT/$PBS_JOBID
TMPDIR=/home/coccia/Codes/workshop_emanuele/WaveT/WaveT/tests/Vacuum_LiCN/field_sse/MPI_tests/$PBS_JOBID
TAPEDIR=$STORAGE/data/$ROOT
mkdir $TMPDIR
#cp $STORAGE/scratch/$ROOT/*inp $TMPDIR
cp /home/coccia/Codes/workshop_emanuele/WaveT/WaveT/tests/Vacuum_LiCN/field_sse/MPI_tests/*inp $TMPDIR 
cd $TMPDIR

EXE=Codes/workshop_emanuele/WaveT/WaveT/bin/WaveT-parallel.x
INP=test_vac_diss.inp
OUT=4out
mpirun -np 4 $HOMEDIR/$EXE < $INP > $OUT
#
# now save appropriate tapes
#
#mkdir $TAPEDIR
#mv *dat $TAPEDIR/
cd ..
#rm -fr $TMPDIR
#ls -l
#
# finished


