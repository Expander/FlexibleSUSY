set terminal epslatex standalone size 15cm,10cm header "\\usepackage{amsmath}"

if (!exists("filename")) filename='higgs-study/data/Alambda/scan_NMSSM_Alambda.dat'

outputfilename=system("echo '".filename."' | sed 's/\\./_/g' ")
outputtexfilename=outputfilename.".tex"

set output outputtexfilename

set title ""
set xlabel '$\tan\beta$'
set xtics format "$%g$"
set ylabel '$m_h^{\text{pole}}$ / GeV'
set ytics format "$%g$"

set key box bottom right width 1 height 0.5

plot "<awk '{ if ($3 == 0) print $0 }' ".filename using 1:2 title '$A^\lambda(M_X)=A_0$', \
     "<awk '{ if ($5 == 0) print $0 }' ".filename using 1:4 title '$A^\lambda(M_X)=10\;\text{TeV}$'     

set output

system "epstopdf ".outputfilename."-inc.eps"
system "pdflatex ".outputtexfilename
