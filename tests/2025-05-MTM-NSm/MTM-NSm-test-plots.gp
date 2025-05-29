# gnuplot file, template generated automatically with
# /Users/gsalam/scripts/gptemplate.py MTM-NSm-test-plots.gp

set term pdfcairo enhanced font 'Helvetica,12pt' lw 1.5 size 10cm,8cm
set  style data lines
set lt 1 pt 4 ps 0.6 lw 1.5 lc rgb '#f00000'
set lt 2 pt 5 ps 0.6 lw 1.5 lc rgb '#00c000'
set lt 3 pt 6 ps 0.6 lw 1.5 lc rgb '#0000e0'
set lt 4 pt 7 ps 0.6 lw 1.5 lc rgb '#ff8000'
set lt 5 pt 8 ps 0.6 lw 1.5 lc rgb '#5070ff'
set key spacing 1.3

set output 'MTM-NSm-test-plots.pdf'

set log x
set xrange [1e-4:]
set xlabel 'x'
set grid

ncol=11
uV="2"
dV="3"
udS="5"
sTot="6"
sV="7"
bTot="10"
g="11"

file2(n)=sprintf("%d",1*n+1*ncol)
uV2=file2(uV)
dV2=file2(dV)
udS2=file2(udS)
sV2=file2(sV)
sTot2=file2(sTot)
bTot2=file2(bTot)
g2=file2(g)

m='<mergeidx.pl '
nnlo_pre =" -f check_nnpdf40_nnlo_preNSm.dat "
nnlo_post=" -f check_nnpdf40_nnlo_postNSm.dat "
nnlo_noev=" -f noev_nnpdf40_nnlo_postNSm.dat "
n3lo_pre =" -f check_nnpdf40_an3lo_preNSm.dat "
n3lo_post=" -f check_nnpdf40_an3lo_postNSm.dat "
n3lo_noev=" -f noev_nnpdf40_an3lo_postNSm.dat "

lhapdf=" LHAPDF.at.Q "
hoppet=" hoppet.at.Q "

nnpdf_N2='NNPDF40\_nnlo\_as\_01180\_mhou'
nnpdf_N3='NNPDF40\_an3lo\_as\_01180'

set label 1 'Q_0 = 4.919 GeV, Q = 4.921 GeV, m_b=4.92 GeV' at graph 0.02,0.95

set title 'NNLO test of b+bbar, '.nnpdf_N2
set ylabel 'x(b + bbar)'
set key bottom
plot m.nnlo_pre .lhapdf u 1:@bTot w l title 'NNLO, LHAPDF', \
     m.nnlo_pre .hoppet u 1:@bTot w l title 'NNLO pre NSm, HOPPET', \
     m.nnlo_post.hoppet u 1:@bTot w l dt 2 title 'NNLO post NSm, HOPPET'

set title 'N3LO test of b+bbar, '.nnpdf_N3
set ylabel 'x(b + bbar)'
set key bottom
plot m.n3lo_pre .lhapdf u 1:@bTot w l title 'N3LO, LHAPDF', \
     m.n3lo_pre .hoppet u 1:@bTot w l title 'N3LO pre NSm, HOPPET', \
     m.n3lo_post.hoppet u 1:@bTot w l dt 2 title 'N3LO post NSm, HOPPET'


set title 'NNLO Δd_V across b threshold, '.nnpdf_N2
set ylabel 'x Δd_V'
plot m.nnlo_noev.nnlo_post.lhapdf u 1:($@dV2-$@dV) w l title 'NNLO, LHAPDF',\
     m.nnlo_noev.nnlo_post.hoppet u 1:($@dV2-$@dV) w l title 'NNLO, HOPPET'

set title 'N3LO Δd_V across b threshold, '.nnpdf_N3
set ylabel 'x Δd_V'
print m.n3lo_noev.n3lo_post.lhapdf
plot m.n3lo_noev.n3lo_post.lhapdf u 1:($@dV2-$@dV) w l title 'N3LO, LHAPDF',\
     m.n3lo_noev.n3lo_pre .hoppet u 1:($@dV2-$@dV) w l title 'N3LO, HOPPET (pre NSm)',\
     m.n3lo_noev.n3lo_post.hoppet u 1:($@dV2-$@dV) w l title 'N3LO, HOPPET (post NSm)'

set title 'N3LO (post-pre)NSm above b threshold, '.nnpdf_N3
set ylabel 'x f(x): post-pre NSm'
plot m.n3lo_pre .n3lo_post.hoppet u 1:($@dV2-$@dV) w l title 'd_V (post-pre)NSm N3LO, HOPPET',\
     m.n3lo_pre .n3lo_post.hoppet u 1:($@uV2-$@uV) w l title 'u_V (post-pre)NSm N3LO, HOPPET',\
     m.n3lo_pre .n3lo_post.hoppet u 1:($@sV2-$@sV) w l title 's_V (post-pre)NSm N3LO, HOPPET',\

set title 'N3LO LHAPDF-HOPPET above b threshold, '.nnpdf_N3
set ylabel 'x f(x): LHAPDF-HOPPET'
plot m.n3lo_pre .hoppet.lhapdf u 1:($@uV2-$@uV) w l title 'u_V N3LO, LHAPDF-HOPPET, pre-NSm',\
     m.n3lo_post.hoppet.lhapdf u 1:($@uV2-$@uV) w l title 'u_V N3LO, LHAPDF-HOPPET, post-NSm N3LO',\

plot m.n3lo_pre .hoppet.lhapdf u 1:($@dV2-$@dV) w l title 'd_V N3LO, LHAPDF-HOPPET, pre-NSm',\
     m.n3lo_post.hoppet.lhapdf u 1:($@dV2-$@dV) w l title 'd_V N3LO, LHAPDF-HOPPET, post-NSm N3LO',\

plot m.n3lo_pre .hoppet.lhapdf u 1:($@udS2-$@udS) w l title 'u_S+d_S N3LO, LHAPDF-HOPPET, pre-NSm',\
     m.n3lo_post.hoppet.lhapdf u 1:($@udS2-$@udS) w l title 'u_S+d_S N3LO, LHAPDF-HOPPET, post-NSm N3LO',\

plot m.n3lo_pre .hoppet.lhapdf u 1:($@uV2+$@dV2+$@udS2-($@uV+$@dV+$@udS)) w l title 'u+ubar+d+dbar N3LO, LHAPDF-HOPPET, pre-NSm',\
     m.n3lo_post.hoppet.lhapdf u 1:($@uV2+$@dV2+$@udS2-($@uV+$@dV+$@udS)) w l title 'u+ubar+d+dbar N3LO, LHAPDF-HOPPET, post-NSm N3LO',\


plot m.n3lo_post.hoppet.lhapdf u 1:(($@g2-$@g)) w l title 'g N3LO, LHAPDF-HOPPET',\


#set title 'N3LO test of d_V, NNPDF40\_an3lo\_as\_01180'
#set ylabel 'x d_V'
#plot m.n3lo_noev.lhapdf u 1:@dV w l title 'N3LO, LHAPDF no ev', \
#     m.n3lo_noev.hoppet u 1:@dV w l title 'N3LO no ev pre NSm, HOPPET', \
#     m.n3lo_post.hoppet u 1:@dV w l dt 2 title 'N3LO post NSm, HOPPET'

# put your code here

set output

