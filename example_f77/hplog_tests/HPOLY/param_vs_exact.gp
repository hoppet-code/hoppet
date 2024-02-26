set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'param_vs_exact.pdf'
set macros

#linetype 1,  linecolor rgb "dark-violet"  linewidth 1.000 dashtype solid pointtype 1 pointsize default
#linetype 2,  linecolor rgb "#009e73"  linewidth 1.000 dashtype solid pointtype 2 pointsize default
#linetype 3,  linecolor rgb "#56b4e9"  linewidth 1.000 dashtype solid pointtype 3 pointsize default
#linetype 4,  linecolor rgb "#e69f00"  linewidth 1.000 dashtype solid pointtype 4 pointsize default
#linetype 5,  linecolor rgb "#f0e442"  linewidth 1.000 dashtype solid pointtype 5 pointsize default
#linetype 6,  linecolor rgb "#0072b2"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
#linetype 7,  linecolor rgb "#e51e10"  linewidth 1.000 dashtype solid pointtype 7 pointsize default
#linetype 8,  linecolor rgb "black" 

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"
ii=1

set key top right

set title 'N3LO'
set ylabel '|C_2^{NS}|'
set xlabel '1-x'

plot 'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:(abs($2))  w l lw 4 title 'XC2NS3PA',\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:(abs($3))  w l lw 4 title 'XC2NS3PE - HPLOG5',\
     'xc2ns3e_vs_xc2ns3p_hpoly.dat' u 1:(abs($3))  w l lw 4 title 'XC2NS3PE - HPOLY'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"
ii=1

set key top right

set title 'N3LO'
set ylabel 'C_2^{NS} param/exact'
set xlabel '1-x'

plot 1 lw 2 lc rgb 'black' not,\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:4  w l lw 4 title 'XC2NS3PA / XC2NS3PE - HPLOG5',\
     'xc2ns3e_vs_xc2ns3p_hpoly.dat' u 1:4  w l lw 4 title 'XC2NS3PA / XC2NS3PE - HPOLY'

set output