set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'
set xlabel 'time (fs)' font "Arial,18"
set output "reaction_field.png"
set ylabel 'Reaction Field (au)' font "Arial,18"
p  "RF_drl_ref/medium_t_1.dat" u ($2*0.024188):5 w l ls 1 lw 1 title "drl",\
   "RF_drl_gen/medium_t_1.dat" u ($2*0.024188):5 w l ls 2 title "gral eps"
set output "local_field.png"
set xrange[500*0.024188:550*0.024188]
set ylabel "Polarization local field (au)" font "Arial,18"
p   "LF_drl_ref/medium_t_1.dat" u ($2*0.024188):5 w l ls 1 title "drl",\
    "LF_drl_gen/medium_t_1.dat" u ($2*0.024188):5 w l ls 2 title "gral eps"
