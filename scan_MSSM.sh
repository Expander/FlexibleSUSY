start=91.1876
stop=100000
n_points=60

./scan.sh --parameter=MS --start=$start --stop=$stop --steps=$n_points --step-size=log --TB=5 --Xt=0 | tee scale_MSSM.dat

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
   12   0                    # force output
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
    0                        # Xt / Ms
"

echo "calculating parametric uncertainty from DeltaMt in FlexibleSUSY/MSSM"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=91~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > scale_MSSMMuBMu_TB-5.dat

{ echo "$slha_templ";
  cat <<EOF
Block EXTPAR
   104  100  # DeltaMt
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=91~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > scale_MSSMMuBMu_TB-5_DeltaMt_high.dat

{ echo "$slha_templ";
  cat <<EOF
Block EXTPAR
   104  -100  # DeltaMt
EOF
} | ./utils/scan-slha.sh \
    --spectrum-generator=models/MSSMMuBMu/run_MSSMMuBMu.x \
    --scan-range=MS[]=91~100000:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    > scale_MSSMMuBMu_TB-5_DeltaMt_low.dat

paste scale_MSSMMuBMu_TB-5.dat \
      scale_MSSMMuBMu_TB-5_DeltaMt_low.dat \
      scale_MSSMMuBMu_TB-5_DeltaMt_high.dat \
      > scale_MSSMMuBMu_TB-5_DeltaMt.dat

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
set style line 14 lt 1 dt 10 lw 1 lc rgb '#FF00FF' pt 3

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

plot [0.091:] [60:140] \
     'scale_MSSM.dat' u (\$1/1000):2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     'scale_MSSM.dat' u (\$1/1000):4 t 'FlexibleSUSY/MSSM' w lines ls 3, \
     'scale_MSSM.dat' u (\$1/1000):5 t 'FlexibleSUSY/HSSUSY' w lines ls 2, \
     'scale_MSSM.dat' u (\$1/1000):6 t 'SOFTSUSY 3.6.2' w lines ls 4, \
     'scale_MSSM.dat' u (\$1/1000):12 t 'SPheno 3.3.8' w lines ls 6, \
     'scale_MSSM.dat' u (\$1/1000):14 t 'FlexibleSUSY/MSSM m_t(M_S)' w points ls 12, \
     'scale_MSSM.dat' u (\$1/1000):15 t 'FlexibleSUSY/MSSM m_t(M_S) SPheno-like' w points ls 13, \
     'scale_MSSM.dat' u (\$1/1000):8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     'scale_MSSM.dat' u (\$1/1000):(\$8-\$9):(\$8+\$9) t 'FeynHiggs uncertainty' w filledcurves ls 9 fs transparent solid 0.3, \
     'scale_MSSM.dat' u (\$1/1000):10 t 'SUSYHD 1.0.2' w lines ls 10, \
     'scale_MSSM.dat' u (\$1/1000):(\$10-\$11):(\$10+\$11) t 'SUSYHD uncertainty' w filledcurves ls 11 fs transparent solid 0.3, \
     'scale_MSSM.dat' u (\$1/1000):16 t 'FlexibleSUSY/MSSM-tower (0L y_t(M_S))' w points ls 14

#     'scale_MSSM.dat' u (\$1/1000):7 t 'FlexibleSUSY/MSSM SPheno-like' w lines ls 5, \
#     'scale_MSSM.dat' u (\$1/1000):13 t 'SPheno/MSSM FS-like' w lines ls 7, \
"

echo "$plot_scale" | gnuplot
