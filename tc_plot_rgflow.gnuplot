set terminal epslatex standalone size 15cm,10cm header "\\usepackage{amsmath}"

outputfilename=system("echo 'tc_gcu' | sed 's/\\./_/g' ")
outputtexfilename=outputfilename.".tex"

set output outputtexfilename

set title ""
set key box top left width -2 height 0.5 spacing 1.5
set xlabel "renormalization scale / GeV"
set logscale x
set format x '$10^{%L}$'
set format y '$%2.1f$'
set mxtics 10

plot [:] \
     filename1 using 1:(column(file1g1column)) title '$g_1^\text{MSSM}$' with lines linewidth 3 linetype 1, \
     filename1 using 1:(column(file1g2column)) title '$g_2^\text{MSSM}$' with lines linewidth 3 linetype 2, \
     filename2 using 1:(column(file2g1column)) title '$g_1^\text{E$_6$SSM}$' with lines linewidth 3 linetype 5, \
     filename2 using 1:(column(file2g2column)) title '$g_2^\text{E$_6$SSM}$' with lines linewidth 3 linetype 4

set output

system "epstopdf ".outputfilename."-inc.eps"
system "pdflatex ".outputtexfilename
