
reset
set st dat li

set log

set xlabel 'evolution time [s]'
set ylabel 'relative accuracy'

set format y "10^{%T}"

set grid

plot '<./cycle2.pl res-pair/dyXX-order6-hires15-dlnlnQXX-du0.4-olnlnQ3.res' u 9:3 w l lw 3 t 'pre-ev, guds, x<0.7 (h15 o:6:4)',\
     '' u 9:6 w l lw 3 t 'pre-ev, gudsc, x<0.9 (h15 o:6:4)',\
     '<./cycle2.pl res-pair/dyXX-order6-hires15-nopreev-dlnlnQXX-du0.4-olnlnQ3.res' u 9:6 w l lw 3 t 'gudsc, x<0.9 (h15 o:6:4)'