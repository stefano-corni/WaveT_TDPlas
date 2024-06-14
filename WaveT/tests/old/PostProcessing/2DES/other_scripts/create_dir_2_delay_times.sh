#!/bin/bash

ICOUNT=50
ISTOP=1000
II=0


C=41.3413745758

while [  $ICOUNT -le $ISTOP ]; do
      JCOUNT=50
      JSTOP=1000
      while [  $JCOUNT -le $JSTOP ]; do
         dir=$ICOUNT'_'$JCOUNT
         mkdir $dir
         cp  inp* ci_* script.sh WaveT-serial.x $dir/
         cd $dir/
         T1=`echo "$ICOUNT*$C" | bc -l` 
         T3=`echo "$JCOUNT*$C" | bc -l` 
         sed -i 's/t1/'"$T1"'/g' inp.inp
         sed -i 's/t3/'"$T3"'/g' inp.inp
         let II=II+1
         sed -i 's/NN/'"$II"'/g' inp.inp 
         cd ../
         let JCOUNT=JCOUNT+50
      done
      let ICOUNT=ICOUNT+50
done

