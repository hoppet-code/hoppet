set terminal pdf enhanced font "Latin Modern Roman,26" size 29cm,23cm
set datafile fortran

set output 'c2_new_fits.pdf'
set macros

f(x)=a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6 + a7*x**7 + a8*(1-x)*log(1-x)
fnf(x)=anf0 + anf1*x + anf2*x**2 + anf3*x**3 + anf4*x**4 + anf5*x**5 + anf6*x**6 + anf7*x**7 
fnf2(x)=anf20 + anf21*x + anf22*x**2 + anf23*x**3 + anf24*x**4 + anf25*x**5 + anf26*x**6 + anf27*x**7 

# Same fit function as 0504242
g(x)=b0 + b1*x + b2*x**2 + b3*x**3 + b4*(1-x)*log(1-x)**2  + log(x)*log(1-x)*(b5+b6*log(x)) + b7*(1-x)*log(1-x)

fit [*:0.97] f(x) 'coefficient_functions_for_fit.dat' i 0 u 1:4 via a0,a1,a2,a3,a4,a5,a6,a7,a8
fit [*:0.97] fnf(x) 'coefficient_functions_for_fit.dat' i 1 u 1:4 via anf0,anf1,anf2,anf3,anf4,anf5,anf6,anf7
fit [*:0.97] fnf2(x) 'coefficient_functions_for_fit.dat' i 2 u 1:4 via anf20,anf21,anf22,anf23,anf24,anf25,anf26,anf27

fit [*:0.97] g(x) 'coefficient_functions_for_fit.dat' i 0 u 1:4 via b0,b1,b2,b3,b4,b5,b6,b7

reset
set mxtics
set mytics
set grid
set log x
#set log y
#set format x "10^{%T}"
#set format y "10^{%T}"
#set yrange [100:*]
ii=1

set key top right

set title 'N3LO'
set ylabel 'C_2^{NS} (exact - expand)'
set xlabel 'x'

plot 'coefficient_functions_for_fit.dat' u 1:4 lw 2 title 'C2NSreg (exact - expand)',\
     f(x) lw 2,\
     g(x) lw 2

reset
set mxtics
set mytics
set grid
set log x
#set log y
#set format x "10^{%T}"
#set format y "10^{%T}"
#set yrange [100:*]
ii=1

set key top right

set title 'N3LO'
set ylabel 'C_2^{NS} (exact - expand)'
set xlabel 'x'

plot 'coefficient_functions_for_fit.dat' u (1.-$1):4 lw 2 title 'C2NSreg (exact - expand)',\
     f(1.-x) lw 2,\
     g(1.-x) lw 2

set output

print a0,'d0 + ', a1, 'd0*y + ', a2, 'd0*y**2 + ', a3, 'd0*y**3 + ', a4, 'd0*y**4 + ', a5, 'd0*y**5 + ', a6, 'd0*y**6 + ', a7, 'd0*y**7'
print ''
print 'nf*(',anf0,'d0 + ', anf1, 'd0*y + ', anf2, 'd0*y**2 + ', anf3, 'd0*y**3 + ', anf4, 'd0*y**4 + ', anf5, 'd0*y**5 + ', anf6, 'd0*y**6 + ', anf7, 'd0*y**7)'
print ''
print 'nf**2*(',anf20,'d0 + ', anf21, 'd0*y + ', anf22, 'd0*y**2 + ', anf23, 'd0*y**3 + ', anf24, 'd0*y**4 + ', anf25, 'd0*y**5 + ', anf26, 'd0*y**6 + ', anf27, 'd0*y**7)'
print ''
print b0,'d0 + ', b1, 'd0*y + ', b2, 'd0*y**2 + ', b3, 'd0*y**3 + ', b4, 'd0*y1*DL1**2 + DL*DL1*(', b5, 'd0 + ', b6, 'd0*DL)'
