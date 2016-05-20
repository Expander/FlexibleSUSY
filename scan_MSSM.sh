start=91.1876
stop=100000
n_points=60

./scan.sh --parameter=MS --start=$start --stop=$stop --steps=$n_points --step-size=log --TB=5 --Xt=0 | tee scale_MSSM.dat

plot_scale="
set terminal pdfcairo
set output 'scale_MSSM.pdf'
set key box bottom right width -4
set logscale x
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

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

plot [0.091:] [60:140] \
     'scale_MSSM.dat' u (\$1/1000):2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     'scale_MSSM.dat' u (\$1/1000):4 t 'FlexibleSUSY/MSSM' w lines ls 3, \
     'scale_MSSM.dat' u (\$1/1000):5 t 'FlexibleSUSY/HSSUSY' w lines ls 2, \
     'scale_MSSM.dat' u (\$1/1000):6 t 'SOFTSUSY 3.6.2' w lines ls 4, \
     'scale_MSSM.dat' u (\$1/1000):7 t 'FlexibleSUSY/MSSM SPheno-like' w lines ls 5, \
     'scale_MSSM.dat' u (\$1/1000):12 t 'SPheno/MSSM' w lines ls 6, \
     'scale_MSSM.dat' u (\$1/1000):13 t 'SPheno/MSSM FS-like' w lines ls 7, \
     'scale_MSSM.dat' u (\$1/1000):14 t 'FlexibleSUSY/MSSM yt(MS)' w points ls 12, \
     'scale_MSSM.dat' u (\$1/1000):15 t 'FlexibleSUSY/MSSM yt(MS) SPheno-like' w points ls 13, \
     'scale_MSSM.dat' u (\$1/1000):8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     'scale_MSSM.dat' u (\$1/1000):(\$8-\$9):(\$8+\$9) t 'FeynHiggs uncertainty' w filledcurves ls 9 fs transparent solid 0.3, \
     'scale_MSSM.dat' u (\$1/1000):10 t 'SUSYHD 1.0.2' w lines ls 10, \
     'scale_MSSM.dat' u (\$1/1000):(\$10-\$11):(\$10+\$11) t 'SUSYHD uncertainty' w filledcurves ls 11 fs transparent solid 0.3
"

echo "$plot_scale" | gnuplot
