set terminal epslatex standalone size 15cm,10cm header "\\usepackage{amsmath}"

outputfilename=system("echo '".filename."' | sed 's/\\./_/g' ")
outputtexfilename=outputfilename.".tex"

set output outputtexfilename

set title ""
set key box bottom right
set xlabel "renormalization scale / GeV"
set logscale x
set format x '$10^{%L}$'
set format y '$%3.2f$'
set mxtics 10

# plot for [i=2:111+1] filename using 1:(column(i)) title columnhead(i)

plot [:] \
     filename using 1:(column(g1column)) title columnhead(g1column) pointtype 4, \
     filename using 1:(column(g2column)) title columnhead(g2column) pointtype 5

set output

system "epstopdf ".outputfilename."-inc.eps"
system "pdflatex ".outputtexfilename
