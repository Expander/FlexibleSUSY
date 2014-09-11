set terminal epslatex standalone size 12cm,8cm header '\renewcommand*{\familydefault}{\sfdefault}\usepackage[cm]{sfmath}\usepackage{amsmath}'

set xlabel "Energieskala / GeV"
set ylabel "Kopplungsst\\\"{a}rke" offset 2
set logscale x
set output "MSSM_gauge_couplings.tex"
set format x "$10^{%L}$"
set key box

if (!exists("filename")) filename='MSSM_rgflow.dat'

set multiplot layout 1,2

set xtics 1,1000,1e18
set size 0.58,1.0
set origin 0.0,0.0

set label 1 "Standardmodell" at 1e6,0.9
set style arrow 1 head filled size screen 0.02,10,45

plot [1:1e17] [0.4:1.2] \
     filename.'.sm' using 1:30 title "$g_1$" with lines lw 3 lt 1, \
     filename.'.sm' using 1:31 title "$g_2$" with lines lw 3 lt 2, \
     filename.'.sm' using 1:32 title "$g_3$" with lines lw 3 lt 3

#set ytics ""
unset key
unset ylabel
set format y ""
set size 0.5,1.0
set origin 0.5,0.0

unset label 1
set label 2 "MSSM" at 1e11,0.9
set label 3 "$\\color{red}{M_X}$" at 2e16,0.28 center
set arrow 1 from 2.e16,0.32 to 2.e16,0.4 as 1 lt rgb "red" lw 2

plot [1:1e17] [0.4:1.2] \
     filename using 1:30 title "$g_1$" with lines lw 3 lt 1, \
     filename using 1:31 title "$g_2$" with lines lw 3 lt 2, \
     filename using 1:32 title "$g_3$" with lines lw 3 lt 3

unset multiplot

set output
system 'pdflatex MSSM_gauge_couplings.tex'
