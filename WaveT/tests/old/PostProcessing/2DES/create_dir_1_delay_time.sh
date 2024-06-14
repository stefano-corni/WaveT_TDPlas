#!/bin/bash

ICOUNT=1
ISTOP=100
II=0
F=4.
START=400.

C=41.3413745758

while [  $ICOUNT -le $ISTOP ]; do
         dir=$ICOUNT'_0'
         mkdir $dir
         cp  inp* ci_* script.sh WaveT-serial.x $dir/
         cd $dir/
         T1=`echo "$START*$C+$ICOUNT*$C/$F" | bc -l` 
         sed -i 's/t1/'"$T1"'/g' inp.inp
         let II=II+1
         sed -i 's/NN/'"$II"'/g' inp.inp 
         cd ../
         let ICOUNT=ICOUNT+1
done

