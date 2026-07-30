#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_mu=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1{a[FNR-1]=$3;next} 
	    FNR!=NR&&FNR>1{d2+=($3-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/mu_t_1.dat mu_t_1.dat)

diff_medium=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1{a[FNR-1]=$3;next} 
	    FNR!=NR&&FNR>1{d2+=($3-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/medium_t_1.dat medium_t_1.dat)


awk 'BEGIN{
      if     ('$diff_mu'< '$threshold0' && '$diff_medium'< '$threshold0'){print  1}
      else if('$diff_mu'> '$threshold0'){print -1}
      else                          {print  0}
   }'
