#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_mu=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/mu_t_1.dat mu_t_1.dat)


diff_ci=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/c_t_1.dat c_t_1.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff_mu'< '$threshold0' && '$diff_ci'< '$threshold0'){print  1}
      else if('$diff_mu'> '$threshold0'){print -1}
      else if('$diff_ci' > '$threshold0'){print -1}
      else                          {print  0}
   }'
