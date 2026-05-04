#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_eps=$(awk '
            FNR==NR&&FNR==2{a=$5;next} 
	    FNR!=NR&&FNR==2{b=$5}
	    END{print sqrt((a-b)*(a-b))}' out/eps eps)


diff_lin=$(awk '
            FNR==NR&&FNR==2{a=$5;next} 
	    FNR!=NR&&FNR==2{b=$5}
	    END{print sqrt((a-b)*(a-b))}' out/lin_entropy lin_entropy)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff_lin'< '$threshold0' && '$diff_eps'< '$threshold0'){print  1}
      else if('$diff_lin'> '$threshold0'){print -1}
      else if('$diff_eps' > '$threshold0'){print -1}
      else                          {print  0}
   }'
