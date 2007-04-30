# gnuplot file

#res-grid/dy0.2-order6-hires15-dlnlnQ0.05-du0.4-olnlnQ4-speed.res

reset

file="res-grid/dy0.2-order6-hires15-dlnlnQ0.05-du0.4-olnlnQ4-speed.res"

set yrange [1:1.2e4]
set ylabel 'Q [GeV]' offset 4
set log y
load "xaxis.gp"
myfont="Helvetica,17"
set xtics font myfont
set ytics font myfont
set ylabel font myfont 
set xlabel font myfont offset 0.0,0.5
set size  0.8

set label 4 "{/ZapfDingbats n} {/Symbol e}&,>10^{-4}" at graph 1.0,1.05 right tc ls 4 font myfont
set label 5 "{/ZapfDingbats n} 10^{-5}<&,{/Symbol e}&,<10^{-4}" at graph 0.52,1.05 tc ls 5 font myfont
set label 6 "{/ZapfDingbats n} 10^{-6}<&,{/Symbol e}&,<10^{-5}" at graph 0.20,1.05 tc ls 6 font myfont
set label 7 "{/ZapfDingbats n} {/Symbol e}&,<10^{-6}" at graph 0.0,1.05 tc ls 7 font myfont

set style line 10 lt 1 pt 5 lc rgb "#c0c0c0"
set style line 4 lt 2 pt 5 ps 0.3
set style line 5 lt 3 pt 5 ps 0.3
set style line 6 lt 1 pt 5 ps 0.3
set style line 7 lt 5 pt 5 ps 0.3



plot \
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 4 -minerr -1 -maxerr -1 -protect" u (zeta_of_y($1)):2 w p ls 10 t '',\
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel -2 -minerr -1 -maxerr -1 -protect" u (zeta_of_y($1)):2 w p ls 10 t '',\
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -minerr 1e-4 -protect" u (zeta_of_y($1)):2 w p ls 4 t '',\
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-4 -minerr 1e-5 -protect" u (zeta_of_y($1)):2 w p ls 5 t '',\
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-5 -minerr 1e-6 -protect" u (zeta_of_y($1)):2 w p ls 6 t '',\
     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-6 -protect" u (zeta_of_y($1)):2 w p ls 7 t ''
