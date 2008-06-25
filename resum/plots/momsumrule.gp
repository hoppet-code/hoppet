set term postscript  enhanced eps color 22 dl 3

set out "momsumrule.eps"

set lmargin 6
set size 2.4,2.2
set origin 0,0

set multiplot

set logscale x
set grid
set xrange [2:1e5]


set xlabel "1/x"






set title "Difference between NLO and NLO + resummation, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 0.2,1.2
set ylabel "xPgg^{res} - xPgg^{nlo}"
set yrange [-0.50:0.50]
set xrange [1:1e7]
unset logscale y

plot "../kfactres-ggnlo-interp.res" u 1:2 t \
 "Interpolation (order 6)" w l lw 4 , \
     "../splresggnlo.res" u 1:6 t "Grid" 

set title "Difference between NLO and NLO + resummation, alphas=0.2, nf=4"
set size 1.0, 1.0
set origin 1.2,1.2
set ylabel "xPqg^{res} - xPqg^{nlo}"
set yrange [-0.50:0.50]
set xrange [1:1e7]
unset logscale y

plot "../kfactres-ggnlo-interp.res" u 1:3 t \
 "Interpolation (order 6)" w l lw 4 , \
     "../splresqgnlo.res" u 1:6 t "Grid" 

set title "Momentum sum rules for gluons"
set size 1.0, 1.0
set origin 0.2,0.2
set ylabel "MSR for gluons"
set xlabel "1/x_{min}"

plot "../msr-nlo.res" u 2:5 t \
"NLO splitting functions"  w l lw 4 , \
"../msr-nlores.res" u 2:5 t \
"Diff (NLOres - NLO)" w l lw 4

set title "Momentum sum rules for quarks"
set size 1.0, 1.0
set origin 1.2,0.2
set ylabel "MSR for quarks"

plot "../msr-nlo.res" u 2:3 t \
"NLO splitting functions"  w l lw 4 , \
"../msr-nlores.res" u 2:3 t \
"Diff (NLOres - NLO)" w l lw 4



reset