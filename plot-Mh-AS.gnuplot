set terminal pdfcairo
filename = "scan_MSSM_AS.dat"
set output filename.".pdf"
set key box top right
set grid

set style line  1 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line  2 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line  3 lt 1 dt 4 lw 2 lc rgb '#45AD53'
set style line  4 lt 1 dt 3 lw 2 lc rgb '#FFBF00'
set style line  5 lt 1 dt 5 lw 2 lc rgb '#FF00FF'
set style line  6 lt 1 dt 6 lw 2 lc rgb '#00FFFF'
set style line  7 lt 1 dt 7 lw 2 lc rgb '#000000'
set style line  8 lt 1 dt 4 lw 2 lc rgb '#00FF00'
set style line  9 lt 1 dt 1 lw 0 lc rgb '#00FF00'
set style line 10 lt 1 dt 2 lw 2 lc rgb '#9C4C17'
set style line 11 lt 1 dt 1 lw 0 lc rgb '#9C4C17'
set style line 12 lt 1 dt 8 lw 1 lc rgb '#000000' pt 1
set style line 13 lt 1 dt 9 lw 1 lc rgb '#000000' pt 2
set style line 14 lt 1 dt 10 lw 1 lc rgb '#FF00FF' pt 3

set xlabel "{/Symbol a}_s(M_Z)"
set ylabel "M_h / GeV"

set arrow from 0.1184,63 to 0.1184,150 nohead lc rgb 'black' dt 2 lw 2
set label at 0.1184,59 center "0.1184"

plot [:] [:] \
     filename u 1:2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     "< awk '{ if ($4 > 0) print }' ".filename u 1:4 t 'FlexibleSUSY/MSSM' w lines ls 3, \
     filename u 1:5 t 'FlexibleSUSY/HSSUSY' w lines ls 2, \
     filename u 1:12 t 'SPheno 3.3.8' w lines ls 6, \
