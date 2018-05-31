set terminal pngcairo
set output "check.png"
set yrange [0:3]
plot "lambda_values.out" u (0):($1):(4):(0) w vectors nohead, "real_eps.inp" u 1:2 w lp

