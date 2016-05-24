start=-3.5
stop=3.5
n_points=60
TB=5
MS=2000

slha_templ="
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   100                  # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   1                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1                    # Higgs logarithmic resummation
   15   1.000000000e-11      # beta-function zero threshold
   16   0                    # calculate observables (a_muon, ...)
   17   0                    # Mt methog (0 = FS)
   18   0                    # print EFT parameters
   19   0                    # mf matching loop order (0 = 1L, 1 = 0L)
Block SMINPUTS               # Standard Model inputs
    0   173.34               # Q_higgs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.384               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
    4   1                    # SignMu
Block EXTPAR
    0   ${MS}                # Qmatch
Block Ms
    ${MS}                    # SUSY scale
Block TanBeta
    ${TB}                    # tan(Beta) at the SUSY scale
Block Xtt
    0                        # Xt / Ms
"

slha_templ_delta_low="
${slha_templ}
Block EXTPAR
    101  -384    # DeltaLambdaASATAT
    102  -384    # DeltaLambdaATATAT
"

slha_templ_delta_high="
${slha_templ}
Block EXTPAR
    101   554.667    # DeltaLambdaASATAT
    102   554.667    # DeltaLambdaATATAT
"

echo "calculating Mh(Xt)"
./scan.sh --parameter=Xt \
          --start=$start \
          --stop=$stop \
          --steps=$n_points \
          --step-size=linear \
          --TB=${TB} \
          --MS=${MS} \
    | tee xt_TB-${TB}_MS-${MS}.dat


echo "calculating parametric uncertainty from Q in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_uncertainty.sh \
        --scan-range=Xtt[]=${start}:${stop}:${n_points} \
        --output=Xtt[],MASS[25] \
        | tee xt_MSSMtower_TB-${TB}_scale_uncertainty.dat

echo "calculating parametric uncertainty from Q_match in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_Qmatch_uncertainty.sh \
        --scan-range=Xtt[]=${start}:${stop}:${n_points} \
        --output=Xtt[],MASS[25] \
        | tee xt_MSSMtower_TB-${TB}_Qmatch_uncertainty.dat

echo "calculating parametric uncertainty from delta in the tower"

echo "$slha_templ_delta_low" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=models/MSSMtower/run_MSSMtower.x \
        --scan-range=Xtt[]=${start}:${stop}:${n_points} \
        --output=Xtt[],MASS[25] \
        | tee xt_MSSMtower_TB-${TB}_delta_low.dat

echo "$slha_templ_delta_high" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=models/MSSMtower/run_MSSMtower.x \
        --scan-range=Xtt[]=${start}:${stop}:${n_points} \
        --output=Xtt[],MASS[25] \
        | tee xt_MSSMtower_TB-${TB}_delta_high.dat

printf "# %14s %16s\n" \
       "MS" "Mh" \
       > xt_MSSMtower_TB-${TB}_Mh.dat

awk '{ if ($1 !~ "#") { printf "%16s %16s\n", $1, $2 } }' xt_TB-${TB}_MS-${MS}.dat \
    >> xt_MSSMtower_TB-${TB}_Mh.dat

paste xt_MSSMtower_TB-${TB}_Mh.dat \
      xt_MSSMtower_TB-${TB}_scale_uncertainty.dat \
      xt_MSSMtower_TB-${TB}_Qmatch_uncertainty.dat \
      xt_MSSMtower_TB-${TB}_delta_low.dat \
      xt_MSSMtower_TB-${TB}_delta_high.dat \
      > xt_MSSMtower_TB-${TB}.dat.$$

printf "# %14s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n" \
       "MS" "Mh" \
       "MS" "Q uncert." \
       "MS" "Qmatch uncert." \
       "MS" "lambda(2L) low" \
       "MS" "lambda(2L) high" \
       > xt_MSSMtower_TB-${TB}_uncertainties.dat

cat xt_MSSMtower_TB-${TB}.dat.$$ >> xt_MSSMtower_TB-${TB}_uncertainties.dat

    plot_scale="
set terminal pdfcairo enhanced size 5in,4in
set output 'xt_TB-${TB}_MS-${MS}.pdf'
set key box top center width -2 at graph 0.5, graph 1.4 opaque
set border back
set grid
set tmargin 8

set style line 1 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line 2 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line 3 lt 1 dt 4 lw 2 lc rgb '#45AD53'
set style line 4 lt 1 dt 3 lw 2 lc rgb '#FFBF00'
set style line 5 lt 1 dt 5 lw 2 lc rgb '#FF00FF'
set style line 6 lt 1 dt 6 lw 2 lc rgb '#00FFFF'
set style line 7 lt 1 dt 1 lw 1 lc rgb '#000000'
set style line 8 lt 1 dt 4 lw 2 lc rgb '#00FF00'
set style line 9 lt 1 dt 1 lw 0 lc rgb '#00FF00'
set style line 10 lt 1 dt 2 lw 2 lc rgb '#9C4C17'
set style line 11 lt 1 dt 1 lw 0 lc rgb '#9C4C17'
set style line 12 lt 1 dt 8 lw 1 lc rgb '#000000' pt 1
set style line 13 lt 1 dt 9 lw 1 lc rgb '#000000' pt 2
set style line 14 lt 1 dt 10 lw 2 lc rgb '#FF00FF' pt 3

set xlabel 'X_t / M_S'
set ylabel 'M_h / GeV'

min(x,y) = x < y ? x : y
max(x,y) = x < y ? y : x

filename = 'xt_TB-${TB}_MS-${MS}.dat'
filename2 = 'xt_MSSMtower_TB-${TB}_uncertainties.dat'

plot [:] [:] \
     filename u 1:2 t 'FS/MSSM-tower' w lines ls 1, \
     filename u 1:4 t 'FS/MSSM' w lines ls 3, \
     filename u 1:5 t 'FS/HSSUSY' w lines ls 2, \
     filename u 1:6 t 'SOFTSUSY 3.6.2' w lines ls 4, \
     filename u 1:7 t 'FS/MSSM SPheno-like' w lines ls 5, \
     filename u 1:8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     filename u 1:10 t 'SUSYHD 1.0.2' w lines ls 10, \
     filename u 1:12 t 'SPheno 3.3.8' w lines ls 6, \
     filename u 1:14 t 'FS/MSSM m_t(M_S)' w points ls 12, \
     filename u 1:15 t 'FS/MSSM m_t(M_S) SPheno-like' w points ls 13, \
     filename u 1:16 t 'FS/MSSM-tower (0L y_t(M_S))' w lines ls 14, \
     filename2 u 1:(\$2-\$4/2):(\$2+\$4/2) t '{/Symbol D} Q' w filledcurves ls 7 fs transparent solid 0.3, \
     filename2 u 1:(\$2-\$6/2):(\$2+\$6/2) t '{/Symbol D} Q_{match}' w filledcurves ls 5 fs transparent solid 0.3, \
     filename2 u 1:(min(\$8,\$10)):(max(\$8,\$10)) t '{/Symbol D} {/Symbol l}^{(2)}' w filledcurves ls 10 fs transparent solid 0.3

"

echo "$plot_scale" | gnuplot
