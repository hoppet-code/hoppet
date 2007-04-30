# gnuplot file

reset
set sty dat li

set ylabel '{/Symbol e}'
set xlabel 'dlnlnQ'
set grid
set log 
set xrange [0.01:0.1]
set yrange [1e-9:1e-2]

set key spacing 1.5
#set key left Left reverse
set key bottom 
set key width -5
set size square
set format y "10^{%T}"

# plots below are given with no-pre, but actually pre or not 
# makes no different here

plot \
     '< ./cycle.pl res-grid/dy0.025-order6-4grids-nopreev-dlnlnQXX-du0.005-olnlnQ3.res' u 1:5 w l lt 2lw 3 t 'all, x<0.9 [o = 3]',\
     '< ./cycle.pl res-grid/dy0.025-order6-4grids-nopreev-dlnlnQXX-du0.005-olnlnQ4.res' u 1:5 w l lt 3 lw 3 t 'all, x<0.9 [o = 4]',\
     ''   u 1:2 w l lt 1 lw 3 t 'guds, x<0.7 [o = 4]'


#     ''   u 1:2 t 'guds x<0.7 [olnlnQ = 3]'