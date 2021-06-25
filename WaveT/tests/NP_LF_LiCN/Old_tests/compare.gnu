set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'

set xlabel 'time (fs)' font "Arial,18"
set ylabel 'NP dipole moment (au)' font "Arial,18"
#set palette defined
set xrange[120:125]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
set style line 1 lt 1 lc rgb "blue"  
set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 

set output "compareDip-Sphere.png"
p  "DIP-Sphere/medium_t_1.dat"     u ($2*0.024188):4 w l ls 1 title "Test",\
   "DIP-Sphere/out/medium_t_1.dat" u ($2*0.024188):4 w l dt 2 title "Reference"

set output "compareOnsOns.png"
p  "ONS-ONS/medium_t_1.dat"     u ($2*0.024188):4 w l ls 1 title "Test",\
   "ONS-ONS/out/medium_t_1.dat" u ($2*0.024188):4 w l dt 2 title "Reference"

set output "compareOnsPCM.png"
p   "ONS-PCM/medium_t_1.dat"     u ($2*0.024188):4 w l ls 1 title "Test",\
    "ONS-PCM/out/medium_t_1.dat" u ($2*0.024188):4 w l dt 2 title "Reference"

set output "compareIEFPCM.png"
p   "IEF-ONS/medium_t_1.dat"     u ($2*0.024188):4 w l ls 1 title "Test",\
    "IEF-ONS/out/medium_t_1.dat" u ($2*0.024188):4 w l dt 2 title "Reference"

set output "compareIEFOns.png"
p   "IEF-ONS/medium_t_1.dat"     u ($2*0.024188):4 w l ls 1 title "Test",\
    "IEF-ONS/out/medium_t_1.dat" u ($2*0.024188):4 w l dt 2 title "Reference"

