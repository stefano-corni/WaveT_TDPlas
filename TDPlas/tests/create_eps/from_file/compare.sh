#!/bin/bash

#gnuplot plot.gnu

t0=0.000000001

d1=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$2;next} 
         FNR!=NR&&FNR>1{d2+=($2-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/eps.out eps.out)

d2=$(awk 'BEGIN{d2=0}
         FNR==NR&&FNR>1{a[FNR-1]=$2;next} 
         FNR!=NR&&FNR>1{d2+=($2-a[FNR-1])^2}
         END{print d2/(FNR-1)}' out/real_imag_eps.out real_imag_eps.out)

awk 'BEGIN{
      if     ('$d1'<'$t0' &&'$d2'<'$t0'){print  1}
      else if('$d1'>'$t0'){print -1}
      else if('$d2'>'$t0'){print -1}
      else                {print  0}
   }'
