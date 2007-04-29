# gnuplot file

reset
set sty dat li

# use any old file
set log x
set log y
set grid
set ylabel 'Q [GeV]'
set xlabel 'x'
set format x "10^{%T}"
set format y "10^{%T}"
set size square
set yrange [1:1.5e4]

plot 'res/dy0.4-order-5-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res' u (exp(-$1)):2 w l lw 3 t ''
