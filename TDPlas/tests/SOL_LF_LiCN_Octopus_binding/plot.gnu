set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'
set output "dipole_vs_t.png"
set xlabel 'Time (fs)' font "Arial,18"
set ylabel 'Dipole in z (D)' font "Arial,18"
p    "Oct-star/td.general/multipoles"        u 1:($4/0.2081943) w l ls 1 title "Octopus w/ internal PCM",\
     "Oct-TDPlas/td.general/multipoles" u 1:($4/0.2081943) w l ls 2 title "Octopus interfaced w/ TDPlas"
