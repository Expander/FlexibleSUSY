start=-3.5
stop=3.5
n_points=60

for MS in 91.1876 200 300 400 500 1000 2000 10000 100000
do
    ./scan.sh --parameter=Xt \
              --start=$start \
              --stop=$stop \
              --steps=$n_points \
              --step-size=linear \
              --TB=5 \
              --MS=${MS} \
        | tee xt_TB-5_MS-${MS}.dat

    plot_scale="
set terminal pdfcairo
set output 'xt_TB-5_MS-${MS}.pdf'
set key box top center width -2
set grid

set style line 1 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line 2 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line 3 lt 1 dt 4 lw 2 lc rgb '#45AD53'
set style line 4 lt 1 dt 3 lw 2 lc rgb '#FFBF00'
set style line 5 lt 1 dt 5 lw 2 lc rgb '#FF00FF'
set style line 6 lt 1 dt 6 lw 1 lc rgb '#00FFFF'
set style line 7 lt 1 dt 1 lw 1 lc rgb '#000000'
set style line 8 lt 1 dt 4 lw 2 lc rgb '#00FF00'
set style line 9 lt 1 dt 1 lw 0 lc rgb '#00FF00'
set style line 10 lt 1 dt 2 lw 2 lc rgb '#9C4C17'
set style line 11 lt 1 dt 1 lw 0 lc rgb '#9C4C17'

set xlabel 'X_t / M_S'
set ylabel 'M_h / GeV'

filename = 'xt_TB-5_MS-${MS}.dat'

plot [:] [:] \
     filename u 1:2 t 'FS/MSSM-tower' w lines ls 1, \
     filename u 1:4 t 'FS/MSSM' w lines ls 3, \
     filename u 1:5 t 'FS/HSSUSY' w lines ls 2, \
     filename u 1:6 t 'SOFTSUSY 3.6.2' w lines ls 4, \
     filename u 1:7 t 'FS/MSSM SPheno-like' w lines ls 5, \
     filename u 1:8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     filename u 1:10 t 'SUSYHD 1.0.2' w lines ls 10
"

    echo "$plot_scale" | gnuplot
done
