set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'
set output "dipole_vs_t.png"
set xlabel 'Time (fs)' font "Arial,18"
set ylabel 'Dipole in z (D)' font "Arial,18"
p "Oct-TDPlas-Sol-only/1-sphere/td.general/multipoles" u ($2*0.6582119569):(sqrt($4**2+$5**2+$6**2)/0.2081943) w l ls 1 title "Solvent",\
  "Oct-TDPlas-Sol-fake-NP/1-sphere/td.general/multipoles" u ($2*0.6582119569):(sqrt($4**2+$5**2+$6**2)/0.2081943) w l ls 2 title "Solvent + NP"
