
reset
set st dat li

set log

set ylabel 'evolution time [s]' offset 2
set xlabel 'relative accuracy'

set format x "10^{%T}"
set xrange [1e-8:1e-2]
set yrange [1e-3:1e1]

set grid

set key spacing 1.5 width -10

set arrow from 5e-6   ,0.00864 to 1.41e-4,0.00864 nohead
set arrow from 1.41e-4,0.00864 to 1.41e-4,0.11 nohead
set label 1 'dy=0.2, dlnlnQ=0.05' at 4e-6,0.008 right

set arrow from 2e-7   ,0.0385 to 6.4e-6,0.0385 nohead
set arrow from 6.4e-6,0.0385  to 6.4e-6,0.2 nohead
set label 2 'dy=0.1, dlnlnQ=0.025' at 5.0e-6,0.28

plot '<./cycle2.pl res-pair/dyXX-order-6-hires15-dlnlnQXX-du0.4-olnlnQ4.res' u 3:9 w l lw 3 t 'pre-ev, guds, x<0.7 (h15 o:-6:4)',\
     '' u 6:9 w l lw 3 t 'pre-ev, gudsc, x<0.9 (h15 o:-6:4)',\
     '<./cycle2.pl res-pair/dyXX-order-6-hires15-nopreev-dlnlnQXX-du0.4-olnlnQ4.res' u 6:9 w l lw 3 t 'gudsc, x<0.9 (h15 o:-6:4)'

#,\
#     '<./cycle2.pl res-pair/dyXX-order-6-hires15-dlnlnQXX-du0.4-olnlnQ4.res' u 8:6 w l lw 3 t 'gudsc, x<0.9 (h15 o:-6:4)'

#     '<./cycle2.pl res-pair/dyXX-order-5-hires15-dlnlnQXX-du0.4-olnlnQ3.res' u 9:3 w l lw 3 t 'pre-ev, guds, x<0.7 (h15 o:6:4)',\
#     '' u 9:6 w l lw 3 t 'pre-ev, gudsc, x<0.9 (h15 o:6:4)'