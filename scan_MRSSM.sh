n_points=60

parameter_point="
Block MINPAR
    3    5                # TanBeta
Block HMIXIN
    301 -0.01             # LSD
    302 -0.01             # LSU
    303 -0.01             # LTD
    304 -0.01             # LTU
    201  1e3              # MuD
    202  1e3              # MuU
"

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
${parameter_point}
"

slha_templ_fs="
${slha_templ}
Block EXTPAR
    0    1e4              # Ms
"

slha_templ_spheno_like="
${slha_templ}
Block FlexibleSUSY
   17   1                  # Mt method (0 = FS, 1 = SPheno)
"

slha_templ_spheno_1L="
${slha_templ}
Block SPhenoInput   # SPheno specific input 
    1  -1              # error level 
    2   0              # SPA conventions 
    7   1              # Skip 2-loop Higgs corrections 
    8   3              # Method used for two-loop calculation 
    9   1              # Gaugeless limit used at two-loop 
   10   1              # safe-mode used at two-loop
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

slha_templ_spheno_2L_AS_low="
${slha_templ_spheno_2L}
Block SMINPUTS               # Standard Model inputs
    3   0.117800             # alpha_s(MZ) SM MSbar
"

slha_templ_spheno_2L_AS_high="
${slha_templ_spheno_2L}
Block SMINPUTS               # Standard Model inputs
    3   0.119000             # alpha_s(MZ) SM MSbar
"

slha_templ_spheno_2L_Mt_low="
${slha_templ_spheno_2L}
Block SMINPUTS               # Standard Model inputs
    6   172.360              # mtop(pole)
"

slha_templ_spheno_2L_Mt_high="
${slha_templ_spheno_2L}
Block SMINPUTS               # Standard Model inputs
    6   174.320              # mtop(pole)
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

slha_templ_AS_low="
${slha_templ}
Block SMINPUTS               # Standard Model inputs
    3   0.117800             # alpha_s(MZ) SM MSbar
"

slha_templ_AS_high="
${slha_templ}
Block SMINPUTS               # Standard Model inputs
    3   0.119000             # alpha_s(MZ) SM MSbar
"

slha_templ_Mt_low="
${slha_templ}
Block SMINPUTS               # Standard Model inputs
    6   172.360              # mtop(pole)
"

slha_templ_Mt_high="
${slha_templ}
Block SMINPUTS               # Standard Model inputs
    6   174.320              # mtop(pole)
"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5.dat

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMMSUSY/run_MRSSMMSUSY.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMMSUSY_TB-5.dat

echo "$slha_templ_spheno_like" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMMSUSY/run_MRSSMMSUSY.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMMSUSY_TB-5_SPheno-like.dat

echo "$slha_templ_spheno_1L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_1L.dat

echo "$slha_templ_spheno_2L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L.dat

echo "$slha_templ_spheno_1L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2_FlexibleSUSY_like \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_1L_FSlike.dat

echo "$slha_templ_spheno_2L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2_FlexibleSUSY_like \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_FSlike.dat

# calculate parametric uncertainty from alpha_s +- 0.0006
echo "calculating parametric uncertainty from alpha_s in SPheno"

echo "$slha_templ_spheno_2L_AS_low" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_AS_low.dat

echo "$slha_templ_spheno_2L_AS_high" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_AS_high.dat

echo "calculating parametric uncertainty from alpha_s in the tower"

echo "$slha_templ_AS_low" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_AS_low.dat

echo "$slha_templ_AS_high" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_AS_high.dat

# calculate parametric uncertainty from Mt +- 0.98 GeV
echo "calculating parametric uncertainty from Mt in SPheno"

echo "$slha_templ_spheno_2L_Mt_low" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_Mt_low.dat

echo "$slha_templ_spheno_2L_Mt_high" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2 \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_Mt_high.dat

echo "calculating parametric uncertainty from Mt in the tower"

echo "$slha_templ_Mt_low" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_Mt_low.dat

echo "$slha_templ_Mt_high" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_Mt_high.dat

# calculate parametric uncertainty from Q
echo "calculating parametric uncertainty from Q in SPheno"

echo "$slha_templ_spheno_2L" | ./utils/scan-slha.sh \
    --spectrum-generator=./SPhenoMRSSM2_uncertainty.sh \
    --scan-range=MINPAR[1]=91~100000:$n_points \
    --step-size=log \
    --output=MINPAR[1],MASS[25] \
    --type=SPheno \
    > scale_SPhenoMRSSM_TB-5_2L_scale_uncertainty.dat

echo "calculating parametric uncertainty from Q in the tower"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=./MRSSMtower_uncertainty.sh \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_scale_uncertainty.dat

# calculate parametric uncertainty from delta lambda +- 1
echo "calculating parametric uncertainty from delta in the tower"

echo "$slha_templ_delta_low" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_delta_low.dat

echo "$slha_templ_delta_high" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:$n_points \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5_delta_high.dat

paste scale_SPhenoMRSSM_TB-5_2L.dat \
      scale_SPhenoMRSSM_TB-5_2L_AS_low.dat \
      scale_SPhenoMRSSM_TB-5_2L_AS_high.dat \
      scale_SPhenoMRSSM_TB-5_2L_Mt_low.dat \
      scale_SPhenoMRSSM_TB-5_2L_Mt_high.dat \
      scale_SPhenoMRSSM_TB-5_2L_scale_uncertainty.dat \
      > scale_SPhenoMRSSM_TB-5_2L.dat.$$

mv scale_SPhenoMRSSM_TB-5_2L.dat.$$ scale_SPhenoMRSSM_TB-5_2L.dat

paste scale_MRSSMtower_TB-5.dat \
      scale_MRSSMtower_TB-5_AS_low.dat \
      scale_MRSSMtower_TB-5_AS_high.dat \
      scale_MRSSMtower_TB-5_Mt_low.dat \
      scale_MRSSMtower_TB-5_Mt_high.dat \
      scale_MRSSMtower_TB-5_delta_low.dat \
      scale_MRSSMtower_TB-5_delta_high.dat \
      scale_MRSSMtower_TB-5_scale_uncertainty.dat \
      > scale_MRSSMtower_TB-5.dat.$$

mv scale_MRSSMtower_TB-5.dat.$$ scale_MRSSMtower_TB-5.dat

# rm -f scale_SPhenoMRSSM_TB-5_2L_Mt_low.dat scale_SPhenoMRSSM_TB-5_2L_Mt_high.dat
# rm -f scale_SPhenoMRSSM_TB-5_2L_AS_low.dat scale_SPhenoMRSSM_TB-5_2L_AS_high.dat
# rm -f scale_MRSSMtower_TB-5_delta_low.dat scale_MRSSMtower_TB-5_delta_high.dat
# rm -f scale_SPhenoMRSSM_TB-5_2L_scale_uncertainty.dat

plot_scale="
set terminal pdfcairo
set output 'scale_MRSSM.pdf'
set key box bottom right width -2
set logscale x
set grid

set style line 1 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line 2 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line 3 lt 1 dt 4 lw 2 lc rgb '#45AD53'
set style line 4 lt 1 dt 3 lw 2 lc rgb '#FFBF00'
set style line 5 lt 1 dt 5 lw 2 lc rgb '#FF00FF'
set style line 6 lt 1 dt 6 lw 2 lc rgb '#00FFFF'
set style line 7 lt 1 dt 7 lw 2 lc rgb '#000000'

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

min(x,y) = x < y ? x : y
max(x,y) = x < y ? y : x

plot [0.09:] [:] \
     'scale_MRSSMtower_TB-5.dat'             u (\$1/1000):2 t 'FS/MRSSM-tower' w lines ls 1, \
     'scale_MRSSMMSUSY_TB-5.dat'             u (\$1/1000):2 t 'FS/MRSSM 1L' w lines ls 3, \
     'scale_MRSSMMSUSY_TB-5_SPheno-like.dat' u (\$1/1000):2 t 'FS/MRSSM 1L SPheno-like' w lines ls 5, \
     'scale_SPhenoMRSSM_TB-5_1L.dat'         u (\$1/1000):2 t 'SPheno/MRSSM 1L' w lines ls 2, \
     'scale_SPhenoMRSSM_TB-5_2L.dat'         u (\$1/1000):2 t 'SPheno/MRSSM 2L' w lines ls 4, \
     'scale_SPhenoMRSSM_TB-5_1L_FSlike.dat'  u (\$1/1000):2 t 'SPheno/MRSSM 1L FS-like' w lines ls 6, \
     'scale_SPhenoMRSSM_TB-5_2L_FSlike.dat'  u (\$1/1000):2 t 'SPheno/MRSSM 2L FS-like' w lines ls 7, \
     'scale_SPhenoMRSSM_TB-5_2L.dat'         u (\$1/1000):(min(\$4,\$6)):(max(\$4,\$6))     t 'alpha_s uncertainty' w filledcurves ls 4 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_SPhenoMRSSM_TB-5_2L.dat'         u (\$1/1000):(min(\$8,\$10)):(max(\$8,\$10))   t 'M_t uncertainty' w filledcurves ls 5 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_SPhenoMRSSM_TB-5_2L.dat'         u (\$1/1000):(\$2-\$12/2):(\$2+\$12/2)         t 'scale uncertainty' w filledcurves ls 6 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_MRSSMtower_TB-5.dat'             u (\$1/1000):(min(\$4,\$6)):(max(\$4,\$6))     t '' w filledcurves ls 4 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_MRSSMtower_TB-5.dat'             u (\$1/1000):(min(\$8,\$10)):(max(\$8,\$10))   t '' w filledcurves ls 5 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_MRSSMtower_TB-5.dat'             u (\$1/1000):(min(\$12,\$14)):(max(\$12,\$14)) t '{/Symbol D}{/Symbol l}^{(2)} uncertainty' w filledcurves ls 1 dt 1 lw 0 fs transparent solid 0.3, \
     'scale_MRSSMtower_TB-5.dat'             u (\$1/1000):(\$2-\$16/2):(\$2+\$16/2)         t '' w filledcurves ls 6 dt 1 lw 0 fs transparent solid 0.3
"

echo "$plot_scale" | gnuplot
