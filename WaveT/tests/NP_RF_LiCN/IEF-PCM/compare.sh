#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001
threshold1=0.0003

diff0_medium=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/medium_t_1.dat medium_t_1.dat)
#echo $diff0_medium

diff1_medium=$(awk '
            NR==2{ print sqrt(($3-$6)^2)}' medium_t_1.dat)


diff_ci=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/c_t_1.dat c_t_1.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff0_medium'< '$threshold0' && '$diff1_medium'< '$threshold1' && '$diff_ci'< '$threshold0'){print  1}
      else if('$diff0_medium'> '$threshold0'){print -1}
      else if('$diff1_medium'> '$threshold1'){print -1}
      else if('$diff_ci' > '$threshold0'){print -1}
      else                          {print  0}
   }'
