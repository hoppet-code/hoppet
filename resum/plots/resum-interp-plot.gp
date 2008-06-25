set term postscript  enhanced eps color 22 dl 3

set out "resum-splitfun-interp.eps"

set lmargin 6
set size 2.4,2.2
set origin 0,0

set multiplot

set logscale x
set grid
set xrange [2:1e5]


set xlabel "1/x"

set title "NLO resummed GG splitting function, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 0.2,0.0
set ylabel "xPgg"
set yrange [0:0.50]
plot "../splresggnlo-interp.res" u 1:2 t "Interpolation (order 5)" w l lw 3 , \
     "../splresggnlo.res" u 1:3 t "Grid" 


set title "NLO resummed GG splitting function, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 1.3,0.0
set ylabel "Relative accuracy"
set yrange [0.0000001:0.50]
set logscale y
plot "../splresggnlo-interp.res" u 1:3 t "Interpolation (order 6) vs. Interpolation (order 5)" 

set title "Difference between NLO and NLO + resummation, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 0.2,1.2
set ylabel "xPgg^{res} - xPgg^{nlo}"
set yrange [-0.50:0.50]
set xrange [1:1e7]
unset logscale y
#plot "../kfatct-ggnlo-interp.res" u 1:2 t \
# "Interpolation (order 5)" w l lw 3 , \
#     "../splresggnlo.res" u 1:6 t "Grid" 

plot "../kfactres-ggnlo-interp.res" u 1:2 t \
 "Interpolation (order 6)" w l lw 4 , \
     "../splresggnlo.res" u 1:6 t "Grid" 


set title "Difference between NLO and NLO + resummation, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 1.3,1.2
set ylabel "Relative accuracy"
set yrange [0.0000001:0.50]
set logscale y
plot "../kfatct-ggnlo-interp.res" u 1:3 t "Interpolation (order 6) vs. Interpolation (order 5)" 



reset