set terminal pdf enhanced font "Latin Modern Roman,32" size 29cm,23cm
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

datafile='coefficient_functions_for_fit.dat'

NF=1 # This is just a dummy to get the coefficient
# The logarithmic large-x coefficients extracted from Large-c-expansion-from-MVV.nb
# C2 coefficients
C2L1x5=-18.962962962962962
C2L1x4=348.44444444444446 
C2L1x3=-1787.041827321787 
C2L1x2=2319.6558207170438 
C2L1x1=5894.634952596364 

C2nfL1x4= - 7.901234567901234*NF
C2nfL1x3=  146.69958847736626*NF 
C2nfL1x2= - 787.5420539087113*NF 
C2nfL1x1=  1199.6906563815378*NF 

C2nf2L1x3=- 0.7901234567901234*NF**2
C2nf2L1x2=  14.617283950617283*NF**2
C2nf2L1x1= - 65.15652656374833*NF**2

# C3 coefficients
C3L1x5=-18.962962962962962
C3L1x4=329.48148148148147
C3L1x3=-1609.6420585153533
C3L1x2=1581.2080878629254
C3L1x1=6889.893780092072

C3nfL1x4= - 7.901234567901234*NF
C3nfL1x3=   134.0576131687243*NF
C3nfL1x2= - 675.1899801124716*NF
C3nfL1x1=   859.3809487061571*NF

C3nf2L1x3=- 0.7901234567901234*NF**2
C3nf2L1x2= 12.246913580246913*NF**2
C3nf2L1x1=- 50.14418088473598*NF**2

# CL coefficients
CLL1x4=18.962962962962962
CLL1x3=-177.39976880643363
CLL1x2=738.4477328541179
CLL1x1=-995.2588274957093

CLnfL1x3= 12.641975308641975*NF
CLnfL1x2=-112.35207379623976*NF
CLnfL1x1=  340.3097076753813*NF

CLnf2L1x2= 2.3703703703703702*NF**2
CLnf2L1x1=-15.012345679012345*NF**2


reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=0

set key at 1e-3,0.65

set title 'C2 N3LO nf independent piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):((C2L1x5*log(1.-$1)**5)/$4) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):((C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4)/$4) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):((C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2+C2L1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=0

set key at 1e-3,1e-4

set title 'C2 N3LO nf independent piece'
set ylabel '|1 - ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C2L1x5*log(1.-$1)**5)/$4)) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4)/$4)) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2+C2L1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=0

set key at 1e-3,1e1

set title 'C2 N3LO nf independent piece'
set ylabel '|Difference wrt C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C2L1x5*log(1.-$1)**5)+$4)) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4)+$4)) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2L1x5*log(1.-$1)**5+C2L1x4*log(1.-$1)**4+C2L1x3*log(1.-$1)**3+C2L1x2*log(1.-$1)**2+C2L1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=1

set key at 1e-3,0.6

set title 'C2 N3LO nf-piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):((C2nfL1x4*log(1.-$1)**4)/$4) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):((C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2+C2nfL1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=1

set key at 1e-3,1e-4

set title 'C2 N3LO nf-piece'
set ylabel '|1 - ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C2nfL1x4*log(1.-$1)**4)/$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2+C2nfL1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=1

set key at 1e-3,1e0

set title 'C2 N3LO nf-piece'
set ylabel '|Difference wrt C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C2nfL1x4*log(1.-$1)**4)+$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2nfL1x4*log(1.-$1)**4+C2nfL1x3*log(1.-$1)**3+C2nfL1x2*log(1.-$1)**2+C2nfL1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=2

set key at 1e-3,0.6

set title 'C2 N3LO nf^2-piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):((C2nf2L1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2+C2nf2L1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=2

set key at 1e-3,1e-3

set title 'C2 N3LO nf^2-piece'
set ylabel '|1 - ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C2nf2L1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2+C2nf2L1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=2

set key at 1e-3,1e-2

set title 'C2 N3LO nf^2-piece'
set ylabel '|Difference wrt C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title 'hep-ph/0504242 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C2nf2L1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C2nf2L1x3*log(1.-$1)**3+C2nf2L1x2*log(1.-$1)**2+C2nf2L1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

# END of C2

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=3

set key at 1e-3,0.65

set title 'C3 N3LO nf independent piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):((C3L1x5*log(1.-$1)**5)/$4) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):((C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4)/$4) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):((C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2+C3L1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=3

set key at 1e-3,1e-4

set title 'C3 N3LO nf independent piece'
set ylabel '|1 - ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C3L1x5*log(1.-$1)**5)/$4)) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4)/$4)) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2+C3L1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=3

set key at 1e-3,1e1

set title 'C3 N3LO nf independent piece'
set ylabel '|Difference wrt C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C3L1x5*log(1.-$1)**5)+$4)) w l lw 6 dt 2 title 'L_5 ln^5(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4)+$4)) w l lw 6 dt 2 title '+ L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3L1x5*log(1.-$1)**5+C3L1x4*log(1.-$1)**4+C3L1x3*log(1.-$1)**3+C3L1x2*log(1.-$1)**2+C3L1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=4

set key at 1e-3,0.6

set title 'C3 N3LO nf-piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):((C3nfL1x4*log(1.-$1)**4)/$4) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):((C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2+C3nfL1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=4

set key at 1e-3,1e-4

set title 'C3 N3LO nf-piece'
set ylabel '|1 - ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C3nfL1x4*log(1.-$1)**4)/$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2+C3nfL1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=4

set key at 1e-3,1e0

set title 'C3 N3LO nf-piece'
set ylabel '|Difference wrt C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C3nfL1x4*log(1.-$1)**4)+$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3nfL1x4*log(1.-$1)**4+C3nfL1x3*log(1.-$1)**3+C3nfL1x2*log(1.-$1)**2+C3nfL1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=5

set key at 1e-3,0.6

set title 'C3 N3LO nf^2-piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):((C3nf2L1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2+C3nf2L1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=5

set key at 1e-3,1e-3

set title 'C3 N3LO nf^2-piece'
set ylabel '|1 - ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(C3nf2L1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2+C3nf2L1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=5

set key at 1e-3,1e-2

set title 'C3 N3LO nf^2-piece'
set ylabel '|Difference wrt C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title '0812.4168 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(C3nf2L1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(C3nf2L1x3*log(1.-$1)**3+C3nf2L1x2*log(1.-$1)**2+C3nf2L1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


# END of C3

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=6

set key at 1e-3,0.65

set title 'CL N3LO nf independent piece'
set ylabel '|Ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):((CLL1x4*log(1.-$1)**4)/$4) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):((CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2+CLL1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=6

set key at 1e-3,1e-4

set title 'CL N3LO nf independent piece'
set ylabel '|1 - ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(CLL1x4*log(1.-$1)**4)/$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2+CLL1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=6

set key at 1e-3,1e1

set title 'CL N3LO nf independent piece'
set ylabel '|Difference wrt C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(CLL1x4*log(1.-$1)**4)+$4)) w l lw 6 dt 2 title 'L_4 ln^4(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title '+ L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLL1x4*log(1.-$1)**4+CLL1x3*log(1.-$1)**3+CLL1x2*log(1.-$1)**2+CLL1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=7

set key at 1e-3,0.6

set title 'CL N3LO nf-piece'
set ylabel '|Ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):((CLnfL1x3*log(1.-$1)**3)/$4) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):((CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2+CLnfL1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=7

set key at 1e-3,1e-4

set title 'CL N3LO nf-piece'
set ylabel '|1 - ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(CLnfL1x3*log(1.-$1)**3)/$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2+CLnfL1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=7

set key at 1e-3,1e0

set title 'CL N3LO nf-piece'
set ylabel '|Difference wrt C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(CLnfL1x3*log(1.-$1)**3)+$4)) w l lw 6 dt 2 title 'L_3 ln^3(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title '+ L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLnfL1x3*log(1.-$1)**3+CLnfL1x2*log(1.-$1)**2+CLnfL1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.2:1.1]
ii=8

set key at 1e-3,0.6

set title 'CL N3LO nf^2-piece'
set ylabel '|Ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):((CLnf2L1x2*log(1.-$1)**2)/$4) w l lw 6 dt 2 title 'L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):((CLnf2L1x2*log(1.-$1)**2+CLnf2L1x1*log(1.-$1)**1)/$4) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)',\
     1 lw 4 lc rgb 'black' not


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=8

set key at 1e-3,1e-3

set title 'CL N3LO nf^2-piece'
set ylabel '|1 - ratio to C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-(CLnf2L1x2*log(1.-$1)**2)/$4)) w l lw 6 dt 2 title 'L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(1.-(CLnf2L1x2*log(1.-$1)**2+CLnf2L1x1*log(1.-$1)**1)/$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

set yrange [*:*]
ii=8

set key at 1e-3,1e-2

set title 'CL N3LO nf^2-piece'
set ylabel '|Difference wrt C_{L,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs($4-$2))  w l lw 6 title ' hep-ph/0411112 parametrisation',\
     datafile i ii u (1.-$1):(abs(-(CLnf2L1x2*log(1.-$1)**2)+$4)) w l lw 6 dt 2 title 'L_2 ln^2(1-x)',\
     datafile i ii u (1.-$1):(abs(-(CLnf2L1x2*log(1.-$1)**2+CLnf2L1x1*log(1.-$1)**1)+$4)) w l lw 6 dt 2 title '+ L_1 ln^1(1-x)'


# END of CL

# START fit

set fit quiet

transition=0.9
a=(transition/(1-transition))**2

damp(x)=a*(1.-x)**2/(a*(1.-x)**2 + x**2)
dummy(x1,x2,x)=x2*damp(x1)+x*(1.-damp(x1))

reset
set dummy x1, x2, x3
fit dummy(x1,x2,x3) datafile i 0 u 1:2:3:4:(($4/$2)**2) zerrors via a

set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=0

set key at 1e-3,0.65

set title 'C2 N3LO nf independent piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+$3*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

#set yrange [0.9:1.1]
ii=0

set key at 1e-3,0.65

set title 'C2 N3LO nf independent piece'
set ylabel '|1- ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-($2*damp($1)+$3*(1.-damp($1)))/$4))  w l lw 6 title 'new parametrisation'

print 'C2 nf indepent: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))

reset
set dummy x1,x2,x3
fit dummy(x1,x2,x3) datafile i 1 u 1:2:3:4:(($4/$2)**2) zerrors via a

set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=1

set key at 1e-3,0.65

set title 'C2 N3LO nf piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+$3*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"

#set yrange [0.9:1.1]
ii=1

set key at 1e-3,0.65

set title 'C2 N3LO nf piece'
set ylabel '|1- ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):(abs(1.-$2/$4))  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(abs(1.-($2*damp($1)+$3*(1.-damp($1)))/$4))  w l lw 6 title 'new parametrisation'

print 'C2 nf piece: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))

reset
set dummy x1, x2, x3
fit dummy(x1,x2,x3) datafile i 2 u 1:2:($3+8149.1250000000000-8070.2796638309956):4:(($4/$2)**2) zerrors via a

set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=2

set key at 1e-3,0.65

set title 'C2 N3LO nf^2 piece'
set ylabel '|Ratio to C_{2,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+($3+8149.1250000000000-8070.2796638309956)*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

print 'C2 nf**2 piece: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))

reset
set dummy x1, x2, x3
fit dummy(x1,x2,x3) datafile i 3 u 1:2:3:4:(($4/$2)**2) zerrors via a
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=3

set key at 1e-3,0.65

set title 'C3 N3LO nf independent piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+$3*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

print 'C3 nf indepent: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))

reset
set dummy x1, x2, x3
fit dummy(x1,x2,x3) datafile i 4 u 1:2:3:4:(($4/$2)**2) zerrors via a
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=3

set key at 1e-3,0.65

set title 'C3 N3LO nf piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+$3*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

print 'C3 nf piece: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))


reset
set dummy x1, x2, x3
fit dummy(x1,x2,x3) datafile i 0 u 1:2:3:4:(($4/$2)**2) zerrors via a
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"
#set format y "10^{%T}"

set yrange [0.9:1.1]
ii=3

set key at 1e-3,0.65

set title 'C3 N3LO nf^2 piece'
set ylabel '|Ratio to C_{3,NS,reg}|'
set xlabel '1-x'

plot datafile i ii u (1.-$1):($2/$4)  w l lw 6 title 'old parametrisation',\
     datafile i ii u (1.-$1):(($2*damp($1)+$3*(1.-damp($1)))/$4)  w l lw 6 title 'new parametrisation',\
     1 lw 4 lc rgb 'black' not

print 'C3 nf**2 piece: a**2 = ', a, ' equivalent t = ', sqrt(a)/(1.+sqrt(a))

set output