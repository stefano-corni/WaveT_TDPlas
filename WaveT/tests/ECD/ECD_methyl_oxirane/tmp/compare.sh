#!/bin/bash

#gnuplot plot.gnu

threshold0=0.0000001

diff_m=$(awk '
            FNR==NR&&FNR==2{a=$3;next}
            FNR!=NR&&FNR==2{b=$3}
            END{print sqrt((a-b)*(a-b))}' out/mu_t_1.dat mu_t_1.dat)


#echo $diff_ci $threshold0

awk 'BEGIN{
      if     ('$diff_m' < '$threshold0') {print  1}
      else if ('$diff_m' > '$threshold0') {print -1}
      else                               {print  0}
   }'
