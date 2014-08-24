set terminal epslatex standalone size 12cm,8cm header '\renewcommand*{\familydefault}{\sfdefault}\usepackage[cm]{sfmath}\usepackage{amsmath}'

set xlabel "Energieskala / GeV"
set logscale x
set output "MSSM_gauge_couplings.tex"
set format x "$10^{%L}$"
set key box

if (!exists("filename")) filename='MSSM_rgflow.dat'

set multiplot layout 1,2

set xtics 1,1000,1e18
set size 0.58,1.0
set origin 0.0,0.0

set label 1 "Standard Model" at 1e6,1

plot [1:1e17] [0.4:1.3] \
     filename.'.sm' using 1:30 title "$g_1$" with lines lw 3 lt 1, \
     filename.'.sm' using 1:31 title "$g_2$" with lines lw 3 lt 2, \
     filename.'.sm' using 1:32 title "$g_3$" with lines lw 3 lt 3

#set ytics ""
unset key
set format y ""
set size 0.5,1.0
set origin 0.5,0.0

unset label 1
set label 2 "MSSM" at 1e11,1

plot [1:1e17] [0.4:1.3] \
     filename using 1:30 title "$g_1$" with lines lw 3 lt 1, \
     filename using 1:31 title "$g_2$" with lines lw 3 lt 2, \
     filename using 1:32 title "$g_3$" with lines lw 3 lt 3

unset multiplot

set output
system 'pdflatex MSSM_gauge_couplings.tex'
