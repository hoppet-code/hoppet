set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'param_vs_exact_C2NS3.pdf'
set macros

#linetype 1,  linecolor rgb "dark-violet"  linewidth 1.000 dashtype solid pointtype 1 pointsize default
#linetype 2,  linecolor rgb "#009e73"  linewidth 1.000 dashtype solid pointtype 2 pointsize default
#linetype 3,  linecolor rgb "#56b4e9"  linewidth 1.000 dashtype solid pointtype 3 pointsize default
#linetype 4,  linecolor rgb "#e69f00"  linewidth 1.000 dashtype solid pointtype 4 pointsize default
#linetype 5,  linecolor rgb "#f0e442"  linewidth 1.000 dashtype solid pointtype 5 pointsize default
#linetype 6,  linecolor rgb "#0072b2"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
#linetype 7,  linecolor rgb "#e51e10"  linewidth 1.000 dashtype solid pointtype 7 pointsize default
#linetype 8,  linecolor rgb "black" 

NF=5
L1x5=-512/27
L1x4=3136/9 - (640* NF)/81
L1x3=-1787.04 + 146.7 *NF - 0.790123 *NF**2
L1x2=2319.66 - 787.542* NF + 14.6173* NF**2
L1x1=4988.49 + 1199.69 *NF - 65.1565 *NF**2

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"
set yrange [100:*]
ii=1

set key top right

set title 'N3LO'
set ylabel '|C_2^{NS}|'
set xlabel '1-x'

plot 'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:(abs($2))  w l lw 4 title 'XC2NS3A_{param}',\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:(abs($3))  w l lw 4 title 'XC2NS3A_{exact}',\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:(L1x5*log($1)**5+L1x4*log($1)**4+L1x3*log($1)**3+L1x2*log($1)**2+L1x1*log($1)**1)  w l lw 4 title 'Log^5[1-x]',\
     
#     'xc2ns3e_vs_xc2ns3p_hpoly.dat' u 1:(abs($3))  w l lw 4 title 'XC2NS3A_{exact} - HPOLY'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=1

set key bottom right

set title 'N3LO'
set ylabel 'C_2^{NS} {param}/{exact}'
set xlabel '1-x'

plot 1 lw 4 lc rgb 'black' not,\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:4  w l lw 4 title 'XC2NS3A_{param} / XC2NS3A_{exact}',\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:((L1x5*log($1)**5+L1x4*log($1)**4+L1x3*log($1)**3+L1x2*log($1)**2+L1x1*log($1)**1)/$3)  w l lw 4 title 'XC2NS3A_{large-x} / XC2NS3A_{exact}',\
     'xc2ns3e_vs_xc2ns3p_hplog.dat' u 1:((L1x5*log($1)**5+(704./3.-640./81.*NF)*log($1)**4+(-3368.+153.5*NF-64./81.*NF**2)*log($1)**3+(-2978.-828.7*NF+18.21*NF**2)*log($1)**2+(18832.-501.1*NF-19.09*NF**2)*log($1)**1)/$3)  w l lw 4 title 'XC2NS3A_{large-x-MVV} / XC2NS3A_{exact}'#,\
#     'xc2ns3e_vs_xc2ns3p_hpoly.dat' u 1:4  w l lw 4 title 'XC2NS3A_{param} / XC2NS3A_{exact} - HPOLY'

set output