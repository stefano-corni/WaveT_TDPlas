#!/bin/bash

#gnuplot plot.gnu

t0=0.0001

d1=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$2;next} 
         FNR!=NR&&FNR>1{d2+=($2-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)

d2=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$3;next} 
         FNR!=NR&&FNR>1{d2+=($3-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)

d3=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$4;next} 
         FNR!=NR&&FNR>1{d2+=($4-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)

d4=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$5;next} 
         FNR!=NR&&FNR>1{d2+=($5-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)

d5=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$6;next} 
         FNR!=NR&&FNR>1{d2+=($6-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)

d6=$(awk 'BEGIN{d2=0}
          FNR==NR&&FNR>1{a[FNR-1]=$7;next} 
	  FNR!=NR&&FNR>1{d2+=($7-a[FNR-1])^2}
	  END{print d2/(FNR-1)}' out/dipole_freq.dat dipole_freq.dat)


awk 'BEGIN{
      if     ('$d1'<'$t0' &&'$d2'<'$t0' &&'$d3'<'$t0' &&'$d4'<'$t0' &&'$d5'<'$t0' &&'$d6'<'$t0'){print  1}
      else if('$d1'>'$t0'){print -1}
      else if('$d2'>'$t0'){print -1}
      else if('$d3'>'$t0'){print -1}
      else if('$d4'>'$t0'){print -1}
      else if('$d5'>'$t0'){print -1}
      else if('$d6'>'$t0'){print -1}
      else                {print  0}
   }'
