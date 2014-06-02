set terminal epslatex standalone size 15cm,10cm header "\\usepackage{amsmath}"

if (!exists("filename")) filename='higgs-study/data/tc/scan_MSSM_tc.dat'

outputfilename=system("echo '".filename."' | sed 's/\\./_/g' ")
outputtexfilename=outputfilename.".tex"

outputGLUfilename=system("echo '".filename."' | sed 's/\\./_Glu_/g' ")
outputGLUtexfilename=outputGLUfilename.".tex"

outputMXfilename=system("echo '".filename."' | sed 's/\\./_MX_/g' ")
outputMXtexfilename=outputMXfilename.".tex"

set output outputtexfilename

set title ""
set xlabel '$t_c$'
set xtics format "$%g$"
set ylabel '$m_h^{\text{pole}}$ / GeV'
set mytics 2

set key box top left width 0.5 height 0.5

plot [:] [121:133] \
     "<awk '{ if ($1 == 5  && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 5\phantom{0}$' , \
     "<awk '{ if ($1 == 10 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 10$', \
     "<awk '{ if ($1 == 20 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 20$', \
     "<awk '{ if ($1 == 30 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 30$', \
     "<awk '{ if ($1 == 40 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 40$'

set output

system "epstopdf ".outputfilename."-inc.eps"
system "pdflatex ".outputtexfilename

### plot gluino mass ###

set output outputGLUtexfilename

set title ""
set xlabel '$t_c$'
set xtics format "$%g$"
set ylabel '$m_{\tilde{g}}^{\text{pole}}$ / TeV'
set mytics 2
set format y '$%2.1f$'

set key box top left width 0.5 height 0.5

plot [:] [:] \
     "<awk '{ if ($1 == 5  && $4 == 0) print $0 }' ".filename using 2:($8/1000) title '$\tan\beta = 5\phantom{0}$' , \
     "<awk '{ if ($1 == 10 && $4 == 0) print $0 }' ".filename using 2:($8/1000) title '$\tan\beta = 10$', \
     "<awk '{ if ($1 == 20 && $4 == 0) print $0 }' ".filename using 2:($8/1000) title '$\tan\beta = 20$', \
     "<awk '{ if ($1 == 30 && $4 == 0) print $0 }' ".filename using 2:($8/1000) title '$\tan\beta = 30$', \
     "<awk '{ if ($1 == 40 && $4 == 0) print $0 }' ".filename using 2:($8/1000) title '$\tan\beta = 40$'

set output

system "epstopdf ".outputGLUfilename."-inc.eps"
system "pdflatex ".outputGLUtexfilename

### plot MX

set output outputMXtexfilename
unset key
set ylabel '$M_X$ / GeV'
set logscale y 10
set mytics 10
set format y '$10^{%L}$'

plot [:] [1.e10:1.e19] \
     "<awk '{ if ($1 == 5  && $4 == 0) print $0 }' ".filename using 2:5 title '' pointtype 5

set output

system "epstopdf ".outputMXfilename."-inc.eps"
system "pdflatex ".outputMXtexfilename
