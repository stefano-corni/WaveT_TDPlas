set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'

set xlabel 'time (fs)' font "Arial,18"
set ylabel 'Reaction Field (au)' font "Arial,18"
#set palette defined
#set xrange[120:125]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
#set style line 1 lt 1 lc rgb "blue"  
#set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 

# from angstrom to atomic units
R = 10.0 * 1.889725989
epsd = 1.18
eps0 = 35.688
tau_deb = 500.0
mu = -3.7219
tau_ons = tau_deb*(2*epsd+1.0)/(2*eps0+1.0)
# from atomic units to femptoseconds
tau_ons_fs = tau_ons*0.0241888425

f(x) = mu/R**3*(2*epsd-2.0)/(2*epsd+1.0)*(1.0-6*(eps0-epsd)/(2*epsd-2)/(2*eps0+1)*(exp(-x/tau_ons_fs)-1.0))
#fit f(x) './PCM/medium_t_1.dat' u ($2*0.024188):6 via R

set output "compareOns.png"
p  "ONS/medium_t_1.dat"     u ($2*0.024188):5 w l ls 1 title "Test",\
   "ONS/medium_t_1.dat" u ($2*0.024188):6 w l dt 2 title "Reference",\
    f(x) title "Onsager analytic"

set output "comparePCM.png"
p   "PCM/medium_t_1.dat"     u ($2*0.024188):5 w l ls 1 title "Test",\
    "PCM/medium_t_1.dat" u ($2*0.024188):6 w l dt 2 title "Reference",\
    f(x) title "Onsager analytic"

set output "comparePCM-mu.png"
p   "PCM-mu/medium_t_1.dat"     u ($2*0.024188):5 w l ls 1 title "Test",\
    "PCM-mu/medium_t_1.dat" u ($2*0.024188):6 w l dt 2 title "Reference"

