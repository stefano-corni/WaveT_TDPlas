#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001
threshold1=0.000000001

diff0=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1{a[FNR-1]=$4;next} 
	    FNR!=NR&&FNR>1{d2+=($4-a[FNR-1])^2}
	    END{print sqrt(d2/(FNR-1))}' out/e_t_1.dat e_t_1.dat)
#echo $diff0 $threshold0
diff1=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1{a[FNR-1]=$7;next} 
	    FNR!=NR&&FNR>1{d2+=($7-a[FNR-1])^2}
	    END{print sqrt(d2/(FNR-1))}' out/e_t_1.dat e_t_1.dat)
#echo $diff1 $threshold1

awk 'BEGIN{
      if     ('$diff0'< '$threshold0' && '$diff1'< '$threshold1'){print  1}
      else if('$diff0'> '$threshold0'){print -1}
      else if('$diff1'> '$threshold1'){print -1}
      else                          {print  0}
   }' 

