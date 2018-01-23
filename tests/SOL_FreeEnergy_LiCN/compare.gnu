set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'

set xlabel 'time (fs)' font "Arial,18"
set ylabel '{/Symbol D}G_{neq} (au)' font "Arial,18"
#set palette defined
#set xrange[120:125]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
#set style line 1 lt 1 lc rgb "blue"  
#set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 

set output "compareOns.png"
p  "ONS/e_t_1.dat"     u ($2*0.024188):6 w l ls 1 title "Test",\
   "ONS/out/e_t_1.dat" u ($2*0.024188):6 w p pt 7 title "Reference"

set output "comparePCM.png"
p   "PCM/e_t_1.dat"     u ($2*0.024188):6 w l ls 1 title "Test",\
    "PCM/out/e_t_1.dat" u ($2*0.024188):6 w p pt 7 title "Reference"

set output "comparePCM-mu.png"
p   "PCM-mu/e_t_1.dat"     u ($2*0.024188):6 w l ls 1 title "Test",\
    "PCM-mu/out/e_t_1.dat" u ($2*0.024188):6 w p pt 7 title "Reference"

