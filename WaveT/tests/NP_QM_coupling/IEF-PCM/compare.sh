#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001
threshold1=0.00000001

diff0_medium=$(awk '
            FNR==NR&&FNR==9{a=$4;next} 
	    FNR!=NR&&FNR==9{b=$4}
	    END{print sqrt((a-b)*(a-b))}' out/g.mat g.mat)

diff1_medium=$(awk '
NR==9{print sqrt(($4-$5)*($4-$5))}' g.mat)

#echo $diff0_medium $threshold0
#echo $diff1_medium $threshold1

awk 'BEGIN{
      if     ('$diff0_medium'< '$threshold0' && '$diff1_medium'< '$threshold1'){print  1}
      else if('$diff0_medium'> '$threshold0'){print -1}
      else if('$diff1_medium'> '$threshold1'){print -1}
      else                          {print  0}
   }'
