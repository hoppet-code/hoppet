# gnuplot file

reset
load 'xaxis.gp'

set log y
set format y "10^{%T}"
set ylabel '{/Symbol e}'
set yrange [1e-7:1e-2]

set key left Left reverse spacing 1.5 width-5

set grid

set label 1 "dy = 0.2" at graph 0.52,0.73
set label 2 "dlnlnQ = 0.05" at graph 0.52,0.66

plot \
 '<./compare2files_v2 res-grid/dy0.2-order-6-4grids-dlnlnQ0.05-du0.4-olnlnQ4.res res-grid/reference.res -channel 11 -protect' u (zeta_of_y($1)):3 w l t 'all flavours',\
 '<./compare2files_v2 res-grid/dy0.2-order-6-4grids-dlnlnQ0.05-du0.4-olnlnQ4.res res-grid/reference.res -channel 4 -minerr -1 -maxerr -1 -protect' u (zeta_of_y($1)):4 w  l lt 3 lw 2 t 'charm near sign-change'

`gnupr acc-dy0.2-dlnlnQ0.05.eps cld`

set yrange [1e-9:1e-4]
set label 1 "dy = 0.05"
set label 2 "dlnlnQ = 0.0125"
plot \
 '<./compare2files_v2 res-grid/dy0.05-order-6-4grids-dlnlnQ0.0125-du0.4-olnlnQ4.res res-grid/reference.res -channel 11 -protect' u (zeta_of_y($1)):3 w l t 'all flavours',\
 '<./compare2files_v2 res-grid/dy0.05-order-6-4grids-dlnlnQ0.0125-du0.4-olnlnQ4.res res-grid/reference.res -channel 4 -minerr -1 -maxerr -1 -protect' u (zeta_of_y($1)):4 w  l lt 3 lw 2 t 'charm near sign-change'

`gnupr acc-dy0.05-dlnlnQ0.0125.eps cld`


