set xlabel 'w1'
set ylabel 'w3'
set view 0,90
set size ratio -1
unset surface
set style data pm3d
set style function pm3d
set pm3d map
set pm3d at b 
set palette model RGB
set palette defined  (-0.001 "white", 1.0 "orange",100. "red",50000. "blue")
#set palette defined (0 "white", 1 "black")
set xrange [0.15:0.55]
set yrange [0.15:0.55]
set zrange [-0.001:50000.00]
sp '2d_spectrum.dat' u 1:2:4
