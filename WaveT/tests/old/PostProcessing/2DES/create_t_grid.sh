#!/bin/bash

ICOUNT=1
ISTOP=100   #number of trajectories
d2=0
dtot=31001  #number of steps
tmid=4000.  #enveloppe center of the first pulse
start=400.  #starting point of the first pulse
dt1=4.0
F=4.

C=41.3413745758

while [  $ICOUNT -le $ISTOP ]; do
         tmp1=`echo "$start*$C+$ICOUNT*$C/$F" | bc ` 
      
         d1=`echo "$tmp1/$dt1" | bc `
         dmid=`echo "$tmid/$dt1" | bc `
         d3=`echo "$dtot - $d1 - $d2 - $dmid" | bc `

         echo $d1 $d3 >> t1_t3.dat

         let ICOUNT=ICOUNT+1
done

