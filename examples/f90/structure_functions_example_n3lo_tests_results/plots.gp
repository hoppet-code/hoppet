set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'plots.pdf'
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
set format x "10^{%T}"
ii=1


set key top left
set title 'F2Z, Q = 100 GeV'

set ylabel 'Ratio to N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 1.0'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-2.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($10+$11+$12+$13)/($23+$24+$25+$26)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 2.0, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-0.5-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($10+$11+$12+$13)/($23+$24+$25+$26)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 0.5, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($10+$11+$12+$13)/($23+$24+$25+$26)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($10+$11+$12+$13)/($23+$24+$25+$26)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 0.5'


set key top right
set title 'F2Z, Q = 100 GeV'

set ylabel 'Ratio to N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 1.0'
set xlabel 'x'

plot 0 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat nnlo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:((($10+$11+$12+$13)-($23+$24+$25))/($10+$11+$12+$13)) lw 2 title 'N^3LO - NNLO, NNLO_{ev}',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-1.0.dat' i ii u 1:((($10+$11+$12+$13)-($23+$24+$25+$26))/($10+$11+$12+$13)) lw 2 title 'N^3LO, NNLO_{ev} - NLO_{ev}'


set key top left
set title 'F2W^+, Q = 100 GeV'

set ylabel 'Ratio to N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 1.0'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-2.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($2+$4+$6+$8)/($15+$17+$19+$21)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 2.0, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-0.5-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($2+$4+$6+$8)/($15+$17+$19+$21)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 0.5, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($2+$4+$6+$8)/($15+$17+$19+$21)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:(($2+$4+$6+$8)/($15+$17+$19+$21)) lw 2 title 'N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 0.5'


set key top right
set title 'F2W^+, Q = 100 GeV'

set ylabel 'Ratio to N^3LO, NNLO_{ev}, x_R = 1.0, x_{R,ev} = 1.0'
set xlabel 'x'

plot 0 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat nnlo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i ii u 1:((($2+$4+$6+$8)-($15+$17+$19))/($2+$4+$6+$8)) lw 2 title 'N^3LO - NNLO, NNLO_{ev}',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-1.0.dat' i ii u 1:((($2+$4+$6+$8)-($15+$17+$19+$21))/($2+$4+$6+$8)) lw 2 title 'N^3LO, NNLO_{ev} - NLO_{ev}'


set key top right
set title 'gluon, Q = 100 GeV'

set ylabel 'Ratio to central NNLO'
set xlabel 'x'
set yrange [0.94:1.14]

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'red' title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'red' title 'N3LO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'red' title 'N3LO_{ev}, x_{R,ev} = 0.5',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'black' title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'black' title 'NNLO_{ev}, x_{R,ev} = 0.5',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'blue' title 'NLO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'blue' title 'NLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($8/$22) lw 2 lc rgb 'blue' title 'NLO_{ev}, x_{R,ev} = 0.5'



set key top right
set title 'u, Q = 100 GeV'

set ylabel 'Ratio to central'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($10/$24) lw 2 title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($10/$24) lw 2 title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($10/$24) lw 2 title 'NNLO_{ev}, x_{R,ev} = 0.5'


set key top right
set title 'd, Q = 100 GeV'

set ylabel 'Ratio to central'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($9/$23) lw 2 title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($9/$23) lw 2 title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($9/$23) lw 2 title 'NNLO_{ev}, x_{R,ev} = 0.5'

set key top right
set title 'ubar, Q = 100 GeV'

set ylabel 'Ratio to central'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($6/$20) lw 2 title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($6/$20) lw 2 title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($6/$20) lw 2 title 'NNLO_{ev}, x_{R,ev} = 0.5'


set key top right
set title 'dbar, Q = 100 GeV'

set ylabel 'Ratio to central'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($7/$21) lw 2 title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($7/$21) lw 2 title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($7/$21) lw 2 title 'NNLO_{ev}, x_{R,ev} = 0.5'


set key top right
set title 's, Q = 100 GeV'

set ylabel 'Ratio to central'
set xlabel 'x'

plot 1 lw 2 lc rgb 'black' notitle,\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($11/$25) lw 2 title 'N3LO_{ev}, x_{R,ev} = 1.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($11/$25) lw 2 title 'NNLO_{ev}, x_{R,ev} = 2.0',\
     '<paste n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat' i 3 u 1:($11/$25) lw 2 title 'NNLO_{ev}, x_{R,ev} = 0.5'




set output