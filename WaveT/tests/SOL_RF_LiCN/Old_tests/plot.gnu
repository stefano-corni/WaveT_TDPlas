set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'
set output "reaction_field.png"
#set palette defined
set xlabel 'time (fs)' font "Arial,18"
set ylabel 'Reaction Field (au)' font "Arial,18"
#set xrange[120:125]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
#set style line 1 lt 1 lc rgb "blue"  
#set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 
p  "ONS/medium_t_1.dat"        u ($2*0.024188):5 w l ls 1 title "Ons",\
   "PCM-mu/medium_t_1.dat" u ($2*0.024188):5 w l ls 2 title "PCM-mu",\
   "PCM/medium_t_1.dat"    u ($2*0.024188):5 w l ls 7 dt 2 title "PCM"

