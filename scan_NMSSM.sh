start=91.1876
stop=100000
n_points=60

DIR=${OUTPUT_DIR:-.}

Xt=0
TB=5
Kappa=0.5
Lambda=0.5


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
Block KappaInput
    ${Kappa}
Block LambdaInput
    ${Lambda}
"

echo "running vanilla NMSSM SGs ..."

echo "$slha_templ" | ./utils/scan-slha.sh \
    --spectrum-generator=models/NMSSMtower/run_NMSSMtower.x \
    --scan-range=MS[]=${start}~${stop}:$n_points \
    --step-size=log \
    --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_NMSSMtower_TB-${TB}.dat

echo "calculating uncertainty from Q_match in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./NMSSMtower_Qmatch_uncertainty.sh \
        --scan-range=MS[]=${start}~${stop}:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_NMSSMtower_TB-${TB}_Xt-${Xt}_Qmatch_uncertainty.dat

echo "uncertainty from Q in the tower"

echo "$slha_templ" | \
    ./utils/scan-slha.sh \
        --spectrum-generator=./NMSSMtower_scale_uncertainty.sh \
        --scan-range=MS[]=${start}~${stop}:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${DIR}"/scale_NMSSMtower_TB-${TB}_Xt-${Xt}_scale_uncertainty.dat

echo "run NMSSM-tower with yt(0L)"

{ echo "$slha_templ";
  cat <<EOF
Block FlexibleSUSY
   19   1    # mf tree-level matching
EOF
} | ./utils/scan-slha.sh \
        --spectrum-generator=models/NMSSMtower/run_NMSSMtower.x \
        --scan-range=MS[]=${start}~${stop}:$n_points \
        --step-size=log \
        --output=MS[],MASS[25] \
    | tee "${Dir}"/scale_NMSSMtower_TB-${TB}_yt-0L.dat
