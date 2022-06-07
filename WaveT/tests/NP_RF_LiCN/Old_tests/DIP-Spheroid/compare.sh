#!/bin/bash

#gnuplot plot.gnu

threshold=0.03

diff=$(awk 'NR==2{print $3-$6}' tmp/medium_t_1.dat)
diff1=$(awk 'NR==2{print $3-$6}' out/medium_t_1.dat)

awk 'BEGIN{
      if     ('$diff'< '$threshold'){print  1}
      else if('$diff'> '$threshold'){print -1}
      else if('$diff'!='$diff1')    {print -1}
      else                          {print  0}
   }' 

