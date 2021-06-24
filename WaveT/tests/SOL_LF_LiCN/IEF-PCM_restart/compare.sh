#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff0_medium=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1&&$1>1980{a[FNR-1]=$5;next} 
	    FNR!=NR&&FNR>1&&$1>1980{d2+=($5-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/medium_t_1.dat medium_t_1.dat)
#echo $diff0_medium


diff_ci=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1&&$1>1980{a[FNR-1]=$3;next} 
            FNR!=NR&&FNR>1&&$1>1980{d2+=($3-a[FNR-1])^2}
            END{print d2/(FNR-1)}' out/c_t_1.dat c_t_1.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff0_medium'< '$threshold0' && '$diff_ci'< '$threshold0'){print  1}
      else if('$diff0_medium'> '$threshold0'){print -1}
      else if('$diff_ci' > '$threshold0'){print -1}
      else                          {print  0}
   }'
