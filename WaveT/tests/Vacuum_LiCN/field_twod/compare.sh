#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_mu_0=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/mu_all_0.dat mu_all_0.dat)

diff_mu_1=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
            FNR!=NR&&FNR==2{b=$3}
            END{print sqrt((a-b)*(a-b))}' out/mu_all_1.dat mu_all_1.dat)


awk 'BEGIN{
      if     ('$diff_mu_0'< '$threshold0' && '$diff_mu_1'< '$threshold0'){print  1}
      else if('$diff_mu_0'> '$threshold0'){print -1}
      else if('$diff_mu_1'> '$threshold0'){print -1}
      else                          {print  0}
   }'
