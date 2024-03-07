set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'param_vs_exact_structure_functions.pdf'
set macros

#linetype 1,  linecolor rgb "dark-violet"  linewidth 1.000 dashtype solid pointtype 1 pointsize default
#linetype 2,  linecolor rgb "#009e73"  linewidth 1.000 dashtype solid pointtype 2 pointsize default
#linetype 3,  linecolor rgb "#56b4e9"  linewidth 1.000 dashtype solid pointtype 3 pointsize default
#linetype 4,  linecolor rgb "#e69f00"  linewidth 1.000 dashtype solid pointtype 4 pointsize default
#linetype 5,  linecolor rgb "#f0e442"  linewidth 1.000 dashtype solid pointtype 5 pointsize default
#linetype 6,  linecolor rgb "#0072b2"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
#linetype 7,  linecolor rgb "#e51e10"  linewidth 1.000 dashtype solid pointtype 7 pointsize default
#linetype 8,  linecolor rgb "black" 

exact='exact_structure_functions_full_exact_tiny_1d-10.dat'
param='param_structure_functions.dat'
#exact='exact_structure_functions_paramNS.dat'

exactparam='<paste '.exact.' '.param

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, W^+'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($8)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($8)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($8)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($8)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($8)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($8)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, W^+'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.6:1.4]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($21/$8) lw 3 title 'F1',\
     exactparam i 1 u 1:($21/$8) lw 3 title 'F2',\
     exactparam i 2 u 1:($21/$8) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO+N3LO, W^+'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($15+$17+$19+$21)/($2+$4+$6+$8)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($15+$17+$19+$21)/($2+$4+$6+$8)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($15+$17+$19+$21)/($2+$4+$6+$8)) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, W^-'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($9)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($9)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($9)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($9)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($9)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($9)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, W^-'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.6:1.4]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($22/$9) lw 3 title 'F1',\
     exactparam i 1 u 1:($22/$9) lw 3 title 'F2',\
     exactparam i 2 u 1:($22/$9) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO+N3LO, W^-'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($16+$18+$20+$22)/($3+$5+$7+$9)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($16+$18+$20+$22)/($3+$5+$7+$9)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($16+$18+$20+$22)/($3+$5+$7+$9)) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, Z'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($13)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($13)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($13)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($13)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($13)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($13)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, N3LO, Z'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.6:1.4]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($26/$13) lw 3 title 'F1',\
     exactparam i 1 u 1:($26/$13) lw 3 title 'F2',\
     exactparam i 2 u 1:($26/$13) lw 3 title 'F3'



reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO+N3LO, Z'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($23+$24+$25+$26)/($10+$11+$12+$13)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($23+$24+$25+$26)/($10+$11+$12+$13)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($23+$24+$25+$26)/($10+$11+$12+$13)) lw 3 title 'F3'

############################################
# NNLO plots below
############################################
reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, W^+'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($6)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($6)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($6)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($6)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($6)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($6)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, W^+'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.96:1.04]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($19/$6) lw 3 title 'F1',\
     exactparam i 1 u 1:($19/$6) lw 3 title 'F2',\
     exactparam i 2 u 1:($19/$6) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO, W^+'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($15+$17+$19)/($2+$4+$6)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($15+$17+$19)/($2+$4+$6)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($15+$17+$19)/($2+$4+$6)) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, W^-'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($7)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($7)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($7)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($7)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($7)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($7)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, W^-'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.96:1.04]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($20/$7) lw 3 title 'F1',\
     exactparam i 1 u 1:($20/$7) lw 3 title 'F2',\
     exactparam i 2 u 1:($20/$7) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO, W^-'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($16+$18+$20)/($3+$5+$7)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($16+$18+$20)/($3+$5+$7)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($16+$18+$20)/($3+$5+$7)) lw 3 title 'F3'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, Z'
set ylabel '|F_i|'
set xlabel 'x'

plot exact i 0 u 1:($1*abs($12)) lw 3 title 'xF_1 exact',\
     param i 0 u 1:($1*abs($12)) lw 3 title 'xF_1 param',\
     exact i 1 u 1:(abs($12)) lw 3 title    'F_2 exact',\
     param i 1 u 1:(abs($12)) lw 3 title    'F_2 param',\
     exact i 2 u 1:($1*abs($12)) lw 3 title 'xF_3 exact',\
     param i 2 u 1:($1*abs($12)) lw 3 title 'xF_3 param'

reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, NNLO, Z'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.96:1.04]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:($25/$12) lw 3 title 'F1',\
     exactparam i 1 u 1:($25/$12) lw 3 title 'F2',\
     exactparam i 2 u 1:($25/$12) lw 3 title 'F3'



reset
set mxtics
set mytics
set grid
set log x
!set log y
set format x "10^{%T}"
ii=1

set key top right

set title 'Q = 100 GeV, LO+NLO+NNLO, Z'
set ylabel 'param/exact'
set xlabel 'x'

set yrange [0.9999:1.0001]

plot 1 lw 2 lc rgb 'black' notitle,\
     exactparam i 0 u 1:(($23+$24+$25)/($10+$11+$12)) lw 3 title 'F1',\
     exactparam i 1 u 1:(($23+$24+$25)/($10+$11+$12)) lw 3 title 'F2',\
     exactparam i 2 u 1:(($23+$24+$25)/($10+$11+$12)) lw 3 title 'F3'
set output