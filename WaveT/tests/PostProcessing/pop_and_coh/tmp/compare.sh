#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_pop=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/pop_t_1.dat pop_t_1.dat)


diff_coh=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/coh_t_1.dat coh_t_1.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff_pop'< '$threshold0' && '$diff_coh'< '$threshold0'){print  1}
      else if('$diff_pop'> '$threshold0'){print -1}
      else if('$diff_coh' > '$threshold0'){print -1}
      else                          {print  0}
   }'
