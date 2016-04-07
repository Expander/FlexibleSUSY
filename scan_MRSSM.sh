slha_templ="
Block MODSEL     #
#   12 1000       # output scale
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
    3    5                # TanBeta
Block EXTPAR
    0    1e4              # Ms
Block MSOFTIN
    300  10               # MDB
    302  2000             # MDO
    301  100              # MDW
Block HMIXIN
    301 -0.1              # LSD
    302 -0.1              # LSU
    303 -0.1              # LTD
    304 -0.1              # LTU
    201  250              # MuD
    202  250              # MuU
    1    0                # Mu
    203  0                # BmuD
    204  0                # BmuU
    101  1e4              # Bmu
"

slha_templ_spheno="
${slha_templ}
Block FlexibleSUSY
   17   1                  # Mt method (0 = FS, 1 = SPheno)
"

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMtower/run_MRSSMtower.x \
    --scan-range=EXTPAR[0]=91~100000:60 \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMtower_TB-5.dat

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMMSUSY/run_MRSSMMSUSY.x \
    --scan-range=EXTPAR[0]=91~100000:60 \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMMSUSY_TB-5.dat

echo "$slha_templ_spheno" | ./utils/scan-slha.sh \
    --spectrum-generator=models/MRSSMMSUSY/run_MRSSMMSUSY.x \
    --scan-range=EXTPAR[0]=91~100000:60 \
    --step-size=log \
    --output=EXTPAR[0],MASS[25] \
    > scale_MRSSMMSUSY_TB-5_SPheno.dat

plot_scale="
set terminal pdfcairo
set output 'scale_MRSSM.pdf'
set key box bottom right
set logscale x
set grid

set style line 1 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line 2 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line 3 lt 1 dt 4 lw 2 lc rgb '#45AD53'
set style line 4 lt 1 dt 3 lw 2 lc rgb '#FFBF00'
set style line 5 lt 1 dt 5 lw 2 lc rgb '#FF00FF'

set xlabel 'M_S / TeV'
set ylabel 'M_h / GeV'

plot [:] [:] \
     'scale_MRSSMtower_TB-5.dat' u (\$1/1000):2 t 'MRSSM-tower' w lines ls 1, \
     'scale_MRSSMMSUSY_TB-5.dat' u (\$1/1000):2 t 'MRSSM 1L' w lines ls 3, \
     'scale_MRSSMMSUSY_TB-5_SPheno.dat' u (\$1/1000):2 t 'MRSSM 1L SPheno-like' w lines ls 5
"

echo "$plot_scale" | gnuplot
