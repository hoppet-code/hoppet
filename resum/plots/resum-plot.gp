set term postscript  enhanced eps color 22 dl 3

set out "resum-splitfun.eps"

set lmargin 6
set size 2.2,2.4
set origin 0,0

set multiplot

set logscale x
set grid
set xrange [2:1e6]


set xlabel "1/x"

set title "GG splitting function, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 0.2,0.0
set ylabel "xPgg"
set yrange [0:0.30]
plot "../splresggnlo.res" u 1:4 t "LO" w l lw 4,\
      "../splresggnlo.res" u 1:5 t "NLO" w l lw 4, \
      "../splresggnlo.res" u 1:3 t "NLO resummed" w l lw 4

set title "QQ splitting function, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 0.2,1.1
set yrange [-0.01:0.10]
set ylabel "xPqq"
plot "../splresqqnlo.res" u 1:4 t "LO" w l lw 4,\
      "../splresqqnlo.res" u 1:5 t "NLO" w l lw 4, \
      "../splresqqnlo.res" u 1:3 t "NLO resummed" w l lw 4

set title "QG splitting function, alphas=0.2, nf=4"
set key left
set size 1.0, 1.0
set origin 1.2,1.1
set yrange [-0.01:0.10]
set ylabel "xPqg"
plot "../splresqgnlo.res" u 1:4 t "LO" w l lw 4,\
      "../splresqgnlo.res" u 1:5 t "NLO" w l lw 4, \
      "../splresqgnlo.res" u 1:3 t "NLO resummed" w l lw 4

set title "GQ splitting function, alphas=0.2, nf=4"
set size 1.0, 1.00
set origin 1.2,0.0
set yrange [-0.00:0.10]
set ylabel "xPgq"
set key bottom
plot "../splresgqnlo.res" u 1:4 t "LO" w l lw 4,\
      "../splresgqnlo.res" u 1:5 t "NLO" w l lw 4, \
      "../splresgqnlo.res" u 1:3 t "NLO resummed" w l lw 4
