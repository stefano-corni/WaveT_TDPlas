#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_mol=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/sp_mol_1.dat sp_mol_1.dat)


diff_np=$(awk '
            FNR==NR&&FNR==2{a=$3;next} 
	    FNR!=NR&&FNR==2{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/sp_np_1.dat sp_np_1.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff_mol'< '$threshold0' && '$diff_np'< '$threshold0'){print  1}
      else if('$diff_mol'> '$threshold0'){print -1}
      else if('$diff_np' > '$threshold0'){print -1}
      else                          {print  0}
   }'
