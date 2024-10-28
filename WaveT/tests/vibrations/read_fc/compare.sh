#!/bin/bash

#gnuplot plot.gnu

threshold0=0.0000001

diff_energy=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1&&$1>10{a[FNR-1]=$4;next} 
	    FNR!=NR&&FNR>1&&$1>10{d2+=($4-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/ci_energy_new.inp ci_energy_new.inp)

diff_mut=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1&&$1>10{a[FNR-1]=$5;next} 
            FNR!=NR&&FNR>1&&$1>10{d2+=($5-a[FNR-1])^2}
            END{print d2/(FNR-1)}' out/ci_mut_new.inp ci_mut_new.inp)

#echo $diff_ci $threshold0

awk 'BEGIN{
      if     ('$diff_energy'< '$threshold0' && '$diff_mut' < '$threshold0'){print  1}
      else if('$diff_energy'> '$threshold0'){print -1}
      else if('$diff_mut' > '$threshold0'){print -1}
      else                          {print  0}
   }'
