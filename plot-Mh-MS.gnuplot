set terminal pdfcairo
set output filename.".pdf"
set key box top left
set logscale x
set grid

set style line 1 lt 1 dt 1 lw 2 lc rgb "#FF0000"
set style line 2 lt 1 dt 2 lw 2 lc rgb "#0000FF"
set style line 3 lt 1 dt 4 lw 2 lc rgb "#45AD53"
set style line 4 lt 1 dt 3 lw 2 lc rgb "#FFBF00"
set style line 5 lt 1 dt 5 lw 2 lc rgb "#FF00FF"

set xlabel "M_S / TeV"
set ylabel "M_h / GeV"

plot [:] [:] \
     filename u ($1/1000):2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     "< awk '{ if ($4 > 0) print }' ".filename u ($1/1000):4 t 'FlexibleSUSY/MSSM 2L' w lines ls 3, \
     filename u ($1/1000):5 t 'HSSUSY 2L' w lines ls 2, \
     filename u ($1/1000):6 t 'SOFTSUSY 2L' w lines ls 4, \
     filename u ($1/1000):7 t 'FlexibleSUSY/MSSM 2L SPheno-like' w lines ls 5
