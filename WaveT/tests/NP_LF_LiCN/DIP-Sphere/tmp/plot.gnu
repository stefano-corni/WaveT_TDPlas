set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'
#set palette defined
set xlabel 'time (fs)' font "Arial,18"
set ylabel 'Dipole moment (au)' font "Arial,18"
set xrange[120:125]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
#set style line 1 lt 1 lc rgb "blue"  
#set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 

set output "IEF-PCM.png"
p   "medium_t_1.dat"    u ($2*0.024188):4 w l ls 1 title "sim",\
   "medium_t_1.dat"    u ($2*0.024188):6 w l ls 2 title "anl",

#set output "DIP-Sphere.png"
#p   "DIP-Sphere/medium_t_1.dat"    u ($2*0.024188):4 w l ls 1 title "sim",\
#   "DIP-Sphere/medium_t_1.dat"    u ($2*0.024188):6 w l ls 2 title "anl",

