
reset
set st dat li

set log

set ylabel 'evolution time [s]' offset 2
set xlabel 'relative accuracy, {/Symbol e}'

set format x "10^{%T}"
set xrange [1e-8:1e-2]
set yrange [1e-3:1e1]

set grid

set key spacing 1.5 width -3

tc=0.0088
acc=7.8736940366E-05
set arrow from acc*0.09,tc to acc,tc nohead
set arrow from acc,tc to acc,tc*10 nohead
set label 1 'dy=0.2, dlnlnQ=0.05' at acc*0.07,tc*0.8 right

tc=0.0389
acc=3.8684380095E-06
set arrow from acc*0.08,tc to acc,tc nohead
set arrow from acc,tc to acc,tc*8 nohead
set label 2 'dy=0.1, dlnlnQ=0.025' at acc,tc*10

plot \
     '<./cycle2.pl res-grid/dyXX-order-6-4grids-nopreev-dlnlnQXX-du0.4-olnlnQ4.res' u 6:9 w l lt 2 lw 3 t 't_i , all, x<0.9',\
     '<./cycle2.pl res-grid/dyXX-order-6-4grids-dlnlnQXX-du0.4-olnlnQ4.res' u 6:9 w l lt 1 lw 3 t 't_c , all, x<0.9',\
     '' u 3:9 w l lt 3 lw 3 t 't_c , guds, x<0.7'

#     '<./cycle2.pl res-pair/dyXX-order-6-hires15-nopreev-dlnlnQXX-du0.4-olnlnQ4.res' u 6:9 w l lw 3 t 'gudsc, x<0.9 (h15 o:-6:4)'

#,\
#     '<./cycle2.pl res-pair/dyXX-order-6-hires15-dlnlnQXX-du0.4-olnlnQ4.res' u 8:6 w l lw 3 t 'gudsc, x<0.9 (h15 o:-6:4)'

#     '<./cycle2.pl res-pair/dyXX-order-5-hires15-dlnlnQXX-du0.4-olnlnQ3.res' u 9:3 w l lw 3 t 'pre-ev, guds, x<0.7 (h15 o:6:4)',\
#     '' u 9:6 w l lw 3 t 'pre-ev, gudsc, x<0.9 (h15 o:6:4)'