reset
set xrange [0.00:  1.00]
set yrange [0.00:  1.00]
#set zrange [0.0:  1.0]
set xtics 0.2
set ytics 0.2
set style data lines
set key left top
set xlabel 'x'
set ylabel 'y'
set view 59,359,1,1
set dgrid3d 41, 41, 2
splot "GNUPLOT_RAMPA2_N1681_E3200.txt" title 'DD method'
pause -1 'Proximo Grafico ?'
