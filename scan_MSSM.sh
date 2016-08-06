start=91.1876
stop=100000
n_points=60

DIR=${OUTPUT_DIR:-.}

# for Xt in 0 2.44949 -2.44949
for Xt in 0 2 -2
do

echo "Running with Xt = ${Xt} (saving in ${DIR})"

./scan.sh --parameter=MS --start=$start --stop=$stop --steps=$n_points --step-size=log --TB=5 --Xt="${Xt}" \
    | tee "${DIR}"/scale_MSSM_TB-5_Xt-${Xt}.dat

slha_templ="
Block MODSEL                 # Select model
    6   0                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
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
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   0.000011663787       # G_Fermi
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
Block TanBeta
    5                        # tan(Beta) at the SUSY scale
Block Xtt
    ${Xt}                    # Xt / Ms
"

echo "calculating parametric uncertainty from DeltaMt in FlexibleSUSY/MSSM"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}.dat

{ echo "$slha_templ";
  cat <<EOF
Block EXTPAR
   104  20.4  # DeltaMt
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_DeltaMt_high.dat

{ echo "$slha_templ";
  cat <<EOF
Block EXTPAR
   104  -20.4  # DeltaMt
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_DeltaMt_low.dat

paste "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}.dat \
      "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_DeltaMt_low.dat \
      "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_DeltaMt_high.dat \
      > "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_DeltaMt.dat

echo "uncertainty from Q in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_uncertainty.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
        | tee "${DIR}"/scale_MSSMtower_TB-${TB}_Xt-${Xt}_scale_uncertainty.dat

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_uncertainty_max.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
        | tee "${DIR}"/scale_MSSMtower_TB-${TB}_Xt-${Xt}_scale_uncertainty_max.dat

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_uncertainty_min.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
        | tee "${DIR}"/scale_MSSMtower_TB-${TB}_Xt-${Xt}_scale_uncertainty_min.dat

echo "calculate Q uncertainty"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=./MSSMMuBMu_uncertainty.sh \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_scale_uncertainty.dat.$$

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=./MSSMMuBMu_uncertainty_max.sh \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_scale_uncertainty_max.dat

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=./MSSMMuBMu_uncertainty_min.sh \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_scale_uncertainty_min.dat

paste "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}.dat \
      "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_scale_uncertainty.dat.$$ \
      > "${DIR}"/scale_MSSMMuBMu_TB-5_Xt-${Xt}_scale_uncertainty.dat

echo "calculating uncertainty from Q_match in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_Qmatch_uncertainty.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_Qmatch_uncertainty.dat

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_Qmatch_uncertainty_max.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_Qmatch_uncertainty_max.dat

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./MSSMtower_Qmatch_uncertainty_min.sh \
        --scan-range=MS[]=${start}~100000:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_Qmatch_uncertainty_min.dat

echo "calculate uncertainty in MSSMtower from varying DeltaLambda"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMtower/run_MSSMtower.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=Ms[],MASS[25] \
    > "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}.dat

{ echo "$slha_templ";
cat <<EOF
Block EXTPAR
    101  -314.485   # DeltaLambdaASATAT
    102  -6.04726   # DeltaLambdaATATAT
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMtower/run_MSSMtower.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=Ms[],MASS[25] \
    > "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_DeltaLambda_low.dat

{ echo "$slha_templ";
cat <<EOF
Block EXTPAR
    101   230.518    # DeltaLambdaASATAT
    102   489.358    # DeltaLambdaATATAT
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMtower/run_MSSMtower.x \
    --scan-range=MS[]=${start}~100000:$n_points \
    --step-size=log \
    --output=Ms[],MASS[25] \
    > "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_DeltaLambda_high.dat

paste "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}.dat \
      "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_DeltaLambda_low.dat \
      "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_DeltaLambda_high.dat \
      > "${DIR}"/scale_MSSMtower_TB-5_Xt-${Xt}_DeltaLambda.dat

plot_scale="
set terminal pdfcairo
set output '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.pdf'
set key box bottom right width -2
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
set style line 14 lt 1 dt 10 lw 1 lc rgb '#FF00FF' pt 3

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

plot [0.091:] [60:140] \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):4 t 'FlexibleSUSY/MSSM' w lines ls 3, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):5 t 'FlexibleSUSY/HSSUSY' w lines ls 2, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):6 t 'SOFTSUSY 3.6.2' w lines ls 4, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):12 t 'SPheno 3.3.8' w lines ls 6, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):(\$8-\$9):(\$8+\$9) t '' w filledcurves ls 9 fs transparent solid 0.3, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):10 t 'SUSYHD 1.0.2' w lines ls 10, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):(\$10-\$11):(\$10+\$11) t '' w filledcurves ls 11 fs transparent solid 0.3, \
     '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):17 t 'SuSpect 2.43' w lines ls 5, \
#    '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):14 t 'FlexibleSUSY/MSSM m_t(M_S)' w points ls 12, \
#    '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):15 t 'FlexibleSUSY/MSSM m_t(M_S) SPheno-like' w points ls 13, \
#    '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):16 t 'FlexibleSUSY/MSSM-tower (0L y_t(M_S))' w points ls 14, \
#    '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):7 t 'FlexibleSUSY/MSSM SPheno-like' w lines ls 5, \
#    '${DIR}/scale_MSSM_TB-5_Xt-${Xt}.dat' u (\$1/1000):13 t 'SPheno/MSSM FS-like' w lines ls 7, \
"

echo "$plot_scale" | gnuplot

done
