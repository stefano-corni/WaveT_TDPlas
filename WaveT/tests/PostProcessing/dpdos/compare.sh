#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001

diff_integ=$(awk '
            FNR==NR&&FNR==30{a=$3;next} 
	    FNR!=NR&&FNR==30{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/integ_td-pdos_frag00001.dat integ_td-pdos_frag00001.dat)


diff_dpdos=$(awk '
            FNR==NR&&FNR==30{a=$3;next} 
	    FNR!=NR&&FNR==30{b=$3}
	    END{print sqrt((a-b)*(a-b))}' out/td-pdos_frag00001.dat td-pdos_frag00001.dat)

#echo $diff_ci

awk 'BEGIN{
      if     ('$diff_integ'< '$threshold0' && '$diff_dpdos'< '$threshold0'){print  1}
      else if('$diff_integ'> '$threshold0'){print -1}
      else if('$diff_dpdos' > '$threshold0'){print -1}
      else                          {print  0}
   }'
