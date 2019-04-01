#!/bin/bash

gnuplot plot.gnu

threshold=0.0001

diff=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1{a[FNR-1]=$5;next} 
	    FNR>1{d2+=($5-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/medium_t_1.dat medium_t_1.dat)

awk 'BEGIN{
      if     ('$diff'< '$threshold'){print  1}
      else if('$diff'> '$threshold'){print -1}
      else                          {print  0}
   }' 

