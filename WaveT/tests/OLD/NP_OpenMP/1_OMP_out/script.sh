#!/bin/sh
#PBS -V
#PBS -N OMP 
#PBS -l nodes=1:ppn=1
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


TMPDIR=$STORAGE/scratch/$ROOT/$PBS_JOBID
TAPEDIR=$STORAGE/data/$ROOT
mkdir $TMPDIR
cp $STORAGE/scratch/$ROOT/*inp $TMPDIR
cd $TMPDIR

EXE=Codes/WaveT/WaveT/bin/WaveT-serial.x
INP=omp.inp
OUT=output
$HOMEDIR/$EXE < $INP > $OUT
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


