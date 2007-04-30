# gnuplot file

reset
set sty dat li

#set ylabel 'rel. precision'
set ylabel '{/Symbol e}'
set xlabel 'dy [base]'
set grid
set log 
set xrange [0.032:0.5]
set mxtics
set yrange [1e-9:1e-2]

set key spacing 1.5 
#set key left Left reverse
set key bottom width -5
set size square
set format y "10^{%T}"

set label 1 'guds, x<0.7' at 0.11,1.5e-7 left rotate by 45
set label 2 'all, x<0.9' at 0.05,2e-6 center rotate by 45

# plots below are given with no-pre, but actually pre or not 
# makes no different here

# 6, -6, -5, -5h9
# 
plot '< ./cycle.pl res-grid/dyXX-order-6-4grids-dlnlnQ0.005-du0.005-olnlnQ4.res' u 1:2 w l lt 1 lw 3 t '4 subgrids (1:3:9:27)',\
     '' u 1:5 w l lt 1 lw 3 t '',\
     '< ./cycle.pl res-grid/dyXX-order-6-hires15-dlnlnQ0.005-du0.005-olnlnQ4.res' u 1:2 w l lt 3 lw 3 t '3 subgrids (1:3:15)',\
     '' u 1:5 w l lt 3 lw 3 t ''

#      '< ./cycle.pl res/dyXX-order-6-hires15-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res' u 1:5 w l lt 2 lw 3 t 'o = -6, h15',\
#      '< ./cycle.pl res/dyXX-order-5-hires15-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res' u 1:5 w l lt 3 lw 3 t 'o = -5, h15',\
#      '< ./cycle.pl res/dyXX-order-5-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res' u 1:5 w l lt 4 lw 3 t 'o = -5, h9'
# 
#     ''   u 1:2 t 'guds x<0.7 [olnlnQ = 3]'