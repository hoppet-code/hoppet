set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'nnpdf_vs_msht.pdf'
set macros

#linetype 1,  linecolor rgb "dark-violet"  linewidth 1.000 dashtype solid pointtype 1 pointsize default
#linetype 2,  linecolor rgb "#009e73"  linewidth 1.000 dashtype solid pointtype 2 pointsize default
#linetype 3,  linecolor rgb "#56b4e9"  linewidth 1.000 dashtype solid pointtype 3 pointsize default
#linetype 4,  linecolor rgb "#e69f00"  linewidth 1.000 dashtype solid pointtype 4 pointsize default
#linetype 5,  linecolor rgb "#f0e442"  linewidth 1.000 dashtype solid pointtype 5 pointsize default
#linetype 6,  linecolor rgb "#0072b2"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
#linetype 7,  linecolor rgb "#e51e10"  linewidth 1.000 dashtype solid pointtype 7 pointsize default
#linetype 8,  linecolor rgb "black" 

nnpdf='n3lo-evolution-nnpdf.dat'
hoppet='n3lo-evolution-hoppet.dat'
msht='n3lo-evolution-msht.dat'
nnlo='nnlo-evolution-hoppet.dat'

all='<paste '.nnpdf.' '.msht.' '.hoppet.' '.nnlo

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel 'x u_v'
set xlabel 'x'

plot all u 1:($20-$29) w l lw 3 title 'Hoppet',\
     all u 1:($11-$29) w l lw 3 title 'MSHT',\
     all u 1:($2-$29) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel 'x d_v'
set xlabel 'x'

plot all u 1:($21-$30) w l lw 3 title 'Hoppet',\
     all u 1:($12-$30) w l lw 3 title 'MSHT',\
     all u 1:($3-$30) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
#set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel 'x (dbar - ubar)'
set xlabel 'x'

plot all u 1:($22-$31) w l lw 3 title 'Hoppet',\
     all u 1:($13-$31) w l lw 3 title 'MSHT',\
     all u 1:($4-$31) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel '|x 2(dbar + ubar)|'
set xlabel 'x'

plot all u 1:(abs($23-$32)) w l lw 3 title 'Hoppet',\
     all u 1:(abs($14-$32)) w l lw 3 title 'MSHT',\
     all u 1:(abs($5-$32)) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel '|x (s + sbar)|'
set xlabel 'x'

plot all u 1:(abs($24-$33)) w l lw 3 title 'Hoppet',\
     all u 1:(abs($15-$33)) w l lw 3 title 'MSHT',\
     all u 1:(abs($6-$33)) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel '|x (c + cbar)|'
set xlabel 'x'

plot all u 1:(abs($25-$34)) w l lw 3 title 'Hoppet',\
     all u 1:(abs($16-$34)) w l lw 3 title 'MSHT',\
     all u 1:(abs($7-$34)) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel '|x (b + bbar)|'
set xlabel 'x'

plot all u 1:(abs($26-$35)) w l lw 3 title 'Hoppet',\
     all u 1:(abs($17-$35)) w l lw 3 title 'MSHT',\
     all u 1:(abs($8-$35)) w l lw 3 title 'NNPDF'

reset
set mxtics
set mytics
set grid
set log x
set log y
set format x "10^{%T}"

set key top right

set title 'Q = 100 GeV, abs. diff. wrt NNLO, VFNS'
set ylabel '|x g|'
set xlabel 'x'

plot all u 1:(abs($27-$36)) w l lw 3 title 'Hoppet',\
     all u 1:(abs($18-$36)) w l lw 3 title 'MSHT',\
     all u 1:(abs($9-$36)) w l lw 3 title 'NNPDF'

set output