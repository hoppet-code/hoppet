# gnuplot file to produce grids of accuracy for Hoppet v2
#
# It assumes that ./run-many.sh as been run to produce a variety of results
# in resdir. 
resdir=" ../../build/pp/"
# It also uses compare2files_v2, for which it is necessary to specify the build dir
builddir="../../build/"


#res-grid/dy0.2-order6-hires15-dlnlnQ0.05-du0.4-olnlnQ4-speed.res

reset
set term pdfcairo enhanced font 'Helvetica,13pt' lw 1.5 size 12cm,8cm
set  style data lines
ps=0.4
set lt 1 ps ps pt 5 lw 1.5 lc rgb '#0000c0'
set lt 2 ps ps pt 5 lw 1.5 lc rgb '#0080ff'
set lt 3 ps ps pt 5 lw 1.5 lc rgb '#008000'
set lt 4 ps ps pt 5 lw 1.5 lc rgb '#ff8000'
set lt 5 ps ps pt 5 lw 1.5 lc rgb '#c00000'
set lt 6 ps ps pt 5 lw 1.5 lc rgb '#f08080'
set lt 7 ps ps pt 5 lw 1.5 lc rgb '#000000'
set key spacing 1.3
set object 1 rectangle from graph 0,0 to graph 1,1 behind fillcolor rgb "grey"

set output "hoppet-v2-accuracy-grids.pdf"


set yrange [1:1.2e4]
set ylabel 'Q [GeV]' offset 2
set log y
load "xaxis.gp"
myfont="Helvetica,13"
set xtics font myfont
set ytics font myfont
set ylabel font myfont 
set xlabel font myfont offset 0.0,0.3
set size  0.78,1

s="◼︎"
#set label 1 s." ε&,<10^{-6}"           at graph -0.03,1.05 tc ls 1 font myfont
#set label 2 s." 10^{-6}<&,ε&,<10^{-5}" at graph  0.10,1.05 tc ls 2 font myfont
#set label 3 s." 10^{-5}<&,ε&,<10^{-4}" at graph  0.30,1.05 tc ls 3 font myfont
#set label 4 s." 10^{-4}<&,ε&,<10^{-3}" at graph  0.50,1.05 tc ls 4 font myfont
#set label 5 s." 10^{-3}<&,ε&,<10^{-2}" at graph  0.70,1.05 tc ls 5 font myfont
#set label 6 s." 10^{-2}<&,ε"           at graph  0.90,1.05 tc ls 6 font myfont

xlab=1.03
set label 1 s." ε&,<10^{-6}"             at graph  xlab,0.90 tc ls 1 font myfont
set label 2 s." 10^{-6}<&,ε&,<&,10^{-5}"   at graph  xlab,0.80 tc ls 2 font myfont
set label 3 s." 10^{-5}<&,ε&,<&,10^{-4}"   at graph  xlab,0.70 tc ls 3 font myfont
set label 4 s." 10^{-4}<&,ε&,<&,10^{-3}"   at graph  xlab,0.60 tc ls 4 font myfont
set label 5 s." 10^{-3}<&,ε&,<&,3×10^{-3}" at graph  xlab,0.50 tc ls 5 font myfont
set label 6 s." 3×10^{-3}<&,ε&,<&,10^{-2}" at graph  xlab,0.40 tc ls 6 font myfont
set label 7 s." 10^{-2}<&,ε"             at graph  xlab,0.30 tc ls 7 font myfont


comp="<".builddir."/benchmarking/compare2files_v2 "
#system("../../build/benchmarking/compare2files_v2")

# as a reference use the highest accuracy grid we have
refrun=resdir."nloop3-ref-dy0.02.dat "

# channel 11 means do all flavours
compargs=" -channel 11 -protect "

set label 101 "NNLO" at graph xlab,0.10
set label 102 "Hoppet v2.0.0" at graph xlab,0.03

dyvals="0.25 0.2 0.15 0.1 0.07 0.05"
#dyvals="0.2 0.15"
do for [dy in dyvals] {
  comprun=resdir."nloop3-preev-dy".dy.".dat "
  #set label 100 "dy = ".dy at graph 1.03,0.9
  set title "Accuracy across all flavours, dy = ".dy
  print "Comparing".refrun." and ".comprun
  plot comp.refrun.comprun.compargs.'-maxerr 1e-6'              u (zeta_of_y($1)):2 w p lt 1 t '',\
       comp.refrun.comprun.compargs.'-maxerr 1e-5 -minerr 1e-6' u (zeta_of_y($1)):2 w p lt 2 t '',\
       comp.refrun.comprun.compargs.'-maxerr 1e-4 -minerr 1e-5' u (zeta_of_y($1)):2 w p lt 3 t '',\
       comp.refrun.comprun.compargs.'-maxerr 1e-3 -minerr 1e-4' u (zeta_of_y($1)):2 w p lt 4 t '',\
       comp.refrun.comprun.compargs.'-maxerr 3e-3 -minerr 1e-3' u (zeta_of_y($1)):2 w p lt 5 t '',\
       comp.refrun.comprun.compargs.'-maxerr 1e-2 -minerr 3e-3' u (zeta_of_y($1)):2 w p lt 6 t '',\
       comp.refrun.comprun.compargs.'             -minerr 1e-2' u (zeta_of_y($1)):2 w p lt 7 t '',\
}


do for [dy in dyvals] {
  comprun=resdir."nloop3-preev-oQ3-oY3-dy".dy.".dat "
  set title "Interp order 3, accuracy across all flavours, dy = ".dy
  print "Comparing".refrun." and ".comprun
  replot
}
do for [dy in dyvals] {
  comprun=resdir."nloop3-preev-oQ2-oY2-dy".dy.".dat "
  set title "Interp order 2, accuracy across all flavours, dy = ".dy
  print "Comparing".refrun." and ".comprun
  replot
}
do for [dy in dyvals] {
  comprun=resdir."nloop3-preev-oQ2-oY3-dy".dy.".dat "
  set title "Interp order y(3),Q(2), accuracy across all flavours, dy = ".dy
  print "Comparing".refrun." and ".comprun
  replot
}


#set title "Exact v. parametrised NNLO splitting+thresholds (all flav)"
#refrun=resdir."nloop3-exactspth-nopreev-dy0.1.dat "
#comprun=resdir."nloop3-nopreev-dy0.1.dat "
#replot

#set title "N3LO v. NNLO splitting"
#refrun=resdir."nloop3-ref-dy0.02.dat "
#comprun=resdir."nloop4-ref-dy0.02.dat "
#replot


set title "Exact v. parametrised NNLO splitting (all flav)"
refrun=resdir."nloop3-exactsp-nopreev-dy0.1.dat "
comprun=resdir."nloop3-nopreev-dy0.1.dat "
replot

set title "Exact v. parametrised NNLO thresholds (all flav)"
refrun =resdir."nloop3-exactspth-nopreev-dy0.1.dat "
comprun=resdir."nloop3-exactsp-nopreev-dy0.1.dat "
replot

# separate the exact-v-param by flavour
flavs="-5 -4 -3 -2 -1 0 1 2 3 4 5"
do for [f in flavs] {
  print "Comparing exact v. param for flavour", f
  set title "Exact v. parametrised NNLO splitting+thresholds (flavour ".f.")"
  #compargs=" -channel ".f." -protect "
  compargs=" -channel ".f." "
  refrun =resdir."nloop3-exactspth-nopreev-dy0.1.dat "
  comprun=resdir."nloop3-nopreev-dy0.1.dat "
  print comp.refrun.comprun.compargs.'-maxerr 1e-4 -minerr 1e-5'
  replot
}

#plot \
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 4 -minerr -1 -maxerr -1 -protect" u (zeta_of_y($1)):2 w p ls 10 t '',\
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel -2 -minerr -1 -maxerr -1 -protect" u (zeta_of_y($1)):2 w p ls 10 t '',\
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -minerr 1e-4 -protect" u (zeta_of_y($1)):2 w p ls 4 t '',\
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-4 -minerr 1e-5 -protect" u (zeta_of_y($1)):2 w p ls 5 t '',\
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-5 -minerr 1e-6 -protect" u (zeta_of_y($1)):2 w p ls 6 t '',\
#     "< ./compare2files_v2 res-grid/dy0.2-order-6-4grids-nopreev-dlnlnQ0.005-du0.005-olnlnQ4.res res-grid/reference.res -channel 11 -maxerr 1e-6 -protect" u (zeta_of_y($1)):2 w p ls 7 t ''
