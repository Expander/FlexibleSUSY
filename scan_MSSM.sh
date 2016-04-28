slha_templ="
Block MODSEL
    1 1           # 1/0: High/low scale input
    2 1           # Boundary Condition
    6 1           # Generation Mixing
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
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
   17   0                    # Mt method (0 = FS, 1 = SPheno)
   18   0                    # write EFT output
Block SMINPUTS               # Standard Model inputs
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
Block MINPAR
    1    1e5              # Ms
    2    0                # Xtt
    3    5                # TanBeta
"

slha_templ_spheno_1L="
${slha_templ}
Block SPhenoInput   # SPheno specific input 
    1  -1              # error level 
    2   0              # SPA conventions 
    7   1              # Skip 2-loop Higgs corrections 
    8   3              # Method used for two-loop calculation 
    9   1              # Gaugeless limit used at two-loop 
   10   0              # safe-mode used at two-loop 
   11   0               # calculate branching ratios 
   13   0               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) 
   14   0               # Run couplings to scale of decaying particle 
   12   1.000E-04       # write only branching ratios larger than this value 
   15   1.000E-30       # write only decay if width larger than this value 
   31   -1              # fixed GUT scale (-1: dynamical GUT scale) 
   32   0               # Strict unification 
   34   1.000E-04       # Precision of mass calculation 
   35   40              # Maximal number of iterations
   36   5               # Minimal number of iterations before discarding points
   37   1               # Set Yukawa scheme  
   38   2               # 1- or 2-Loop RGEs 
   50   1               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) 
   51   0               # Write Output in CKM basis 
   52   0               # Write spectrum in case of tachyonic states 
   55   1               # Calculate loop corrected masses 
   57   0               # Calculate low energy constraints 
   65   1               # Solution tadpole equation 
   75   1               # Write WHIZARD files 
   76   1               # Write HiggsBounds file   
   86   0               # Maximal width to be counted as invisible in Higgs decays; -1: only LSP 
  510   0               # Write tree level values for tadpole solutions 
  515   0               # Write parameter values at GUT scale 
  520   0               # Write effective Higgs couplings (HiggsBounds blocks): put 0 to use file with MadGraph! 
  521   0               # Diphoton/Digluon widths including higher order 
  525   0               # Write loop contributions to diphoton decay of Higgs 
  530   0               # Write Blocks for Vevacious 
"

slha_templ_spheno_2L="
${slha_templ_spheno_1L}
Block SPhenoInput   # SPheno specific input 
    7   0              # Skip 2-loop Higgs corrections 
"

start=91.1876
stop=100000
n_points=60

./scan.sh --parameter=MS --start=$start --stop=$stop --steps=$n_points --step-size=log --TB=5 --Xt=0 | tee scale_MSSM.dat

echo "$slha_templ_spheno_2L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMSSM \
    --scan-range=MINPAR[1]=$start~$stop:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    | tee scale_SPhenoMSSM_TB-5_2L.dat

echo "$slha_templ_spheno_2L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMSSM_FlexibleSUSY_like \
    --scan-range=MINPAR[1]=$start~$stop:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    | tee scale_SPhenoMSSM_TB-5_2L_FSlike.dat

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

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

plot [:] [60:140] \
     'scale_MSSM.dat'                      u (\$1/1000):2 t 'FlexibleSUSY/MSSM-tower' w lines ls 1, \
     'scale_MSSM.dat'                      u (\$1/1000):4 t 'FlexibleSUSY/MSSM 2L' w lines ls 3, \
     'scale_MSSM.dat'                      u (\$1/1000):5 t 'HSSUSY 2L' w lines ls 2, \
     'scale_MSSM.dat'                      u (\$1/1000):6 t 'SOFTSUSY 2L' w lines ls 4, \
     'scale_MSSM.dat'                      u (\$1/1000):7 t 'FlexibleSUSY/MSSM 2L SPheno-like' w lines ls 5, \
     'scale_SPhenoMSSM_TB-5_2L.dat'        u (\$1/1000):2 t 'SPheno/MSSM 2L' w lines ls 6, \
     'scale_SPhenoMSSM_TB-5_2L_FSlike.dat' u (\$1/1000):2 t 'SPheno/MSSM 2L FS-like' w lines ls 7, \
     'scale_MSSM.dat'                      u (\$1/1000):8 t 'FeynHiggs 2.11.3' w lines ls 8, \
     'scale_MSSM.dat'                      u (\$1/1000):(\$8-\$9):(\$8+\$9) t 'FeynHiggs uncertainty' w filledcurves ls 9 fs transparent solid 0.3, \
     'scale_MSSM.dat'                      u (\$1/1000):10 t 'SUSYHD 1.0.2' w lines ls 10, \
     'scale_MSSM.dat'                      u (\$1/1000):(\$10-\$11):(\$10+\$11) t 'SUSYHD uncertainty' w filledcurves ls 11 fs transparent solid 0.3
"

echo "$plot_scale" | gnuplot
