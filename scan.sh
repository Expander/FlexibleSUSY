#!/bin/sh

output="MASS-25"
parameter=MS
Xt=0
TB=5
MS=100000
At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
M3factor=0.99999
M3=
AS="1.184000000e-01"
AI="1.279440000e+02"
MT="1.733400000e+02"
MTmethod=0
WRITE_EFT=0
MF_TL_MATCHING=0
GF=0.0000116638
MZ=91.1876
BL=3

dump_fs_slha_input_file=
dump_fs_slha_output_file=
dump_ss_slha_input_file=
dump_ss_slha_output_file=

start=91
stop=1000
steps=10

step_size=linear

sminputs_tmpl="\
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
"

slha_tmpl="\
Block MODSEL                 # Select model
    6   0                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   1000                 # max. iterations (0 = automatic)
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
${sminputs_tmpl}
Block MINPAR                 # Input parameters
    4   1                    # SignMu
Block EXTPAR                 # Input parameters
  100   2                    # LambdaLoopOrder (HSSUSY)
Block Ms
    ${MS}                    # SUSY scale
Block TanBeta
    ${TB}                    # tan(Beta) at the SUSY scale
Block Xtt
    ${Xt}                    # Xt / Ms
"

# prints SLHA block
print_slha_block_awk='
BEGIN {
   is_block = 0;
   if (block == "") {
      print "Error: block name not defined";
      print "   Please define the block name with -v block=<block-name>";
      exit 1
   }
}
{
   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$";
   not_pattern = "^block[[:blank:]]*.*$";

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   };

   if (is_block)
      print $0
}
'

# prints block entry
# expects block entry keys in the form x or x:y or x:y:z etc.
print_block_entry_awk='
{
  len = split(keys,k,":");

  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(len + 1)
}
'

run_sg() {
    local SG="$1"
    local MS2=$(echo "scale=5; ${MS}^2" | bc)
    local At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    local slha_output=
    local block=
    local value=
    local slha_input=
    local output_block=$(echo "${output}" | cut -d'-' -f1)
    local output_entry=$(echo "${output}" | cut -d'-' -f2)

    slha_input=$(
    { echo "$slha_tmpl" ; \
      cat <<EOF
Block MODSEL                 # Select model
   12    ${MS}
Block FlexibleSUSY
    6   ${BL}                # beta-functions loop order
   17   ${MTmethod}          # mt calculation (0 = FlexibleSUSY, 1 = SPheno)
   18   ${WRITE_EFT}         # write full model (0) / EFT (1)
   19   ${MF_TL_MATCHING}    # mf tree-level matching
Block SMINPUTS               # Standard Model inputs
    1   ${AI}                # alpha_em(MZ) SM MSbar
    2   ${GF}                # G_Fermi
    3   ${AS}                # alpha_s(MZ) SM MSbar
    4   ${MZ}                # MZ(pole)
    6   ${MT}                # mtop(pole)
Block TanBeta
    ${TB}                    # tan(Beta) at the SUSY scale
Block Xtt
    ${Xt}                    # Xt / Ms
Block Ms
    ${MS}  # SUSY scale
Block EXTPAR                 # Input parameters
    0   ${MS}                # MSUSY
    1   ${MS}                # M1(MSUSY)
    2   ${MS}                # M2(MSUSY)
    3   ${M3}                # M3(MSUSY)
    4   ${MS}                # Mu(MSUSY)
    5   ${MS}                # mA(MSUSY)
    6   173.34               # MEWSB
    7   ${At}                # At(MSUSY)
   25   ${TB}                # TanBeta(MSUSY)
Block MSQ2IN
  1  1     ${MS2}   # mq2(1,1)
  2  2     ${MS2}   # mq2(2,2)
  3  3     ${MS2}   # mq2(3,3)
Block MSE2IN
  1  1     ${MS2}   # me2(1,1)
  2  2     ${MS2}   # me2(2,2)
  3  3     ${MS2}   # me2(3,3)
Block MSL2IN
  1  1     ${MS2}   # ml2(1,1)
  2  2     ${MS2}   # ml2(2,2)
  3  3     ${MS2}   # ml2(3,3)
Block MSU2IN
  1  1     ${MS2}   # mu2(1,1)
  2  2     ${MS2}   # mu2(2,2)
  3  3     ${MS2}   # mu2(3,3)
Block MSD2IN
  1  1     ${MS2}   # md2(1,1)
  2  2     ${MS2}   # md2(2,2)
  3  3     ${MS2}   # md2(3,3)
EOF
    })

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    [ "x$dump_fs_slha_input_file" != "x" ] && \
        echo "$slha_input" > "$dump_fs_slha_input_file"

    [ "x$dump_fs_slha_output_file" != "x" ] && \
        echo "$slha_output" > "$dump_fs_slha_output_file"

    echo $value
}

run_spheno() {
    local SG="$1"
    local MS2=$(echo "scale=5; ${MS}^2" | bc)
    local At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    local slha_output=
    local block=
    local value=
    local slha_input=
    local output_block=$(echo "${output}" | cut -d'-' -f1)
    local output_entry=$(echo "${output}" | cut -d'-' -f2)

    cat <<EOF > SPheno.in
Block MODSEL                 # Select model
    1 1           # 1/0: High/low scale input
    2 1           # Boundary Condition
   12    ${MS}
Block SMINPUTS               # Standard Model inputs
    1   ${AI}                # alpha_em(MZ) SM MSbar
    2   ${GF}                # G_Fermi
    3   ${AS}                # alpha_s(MZ) SM MSbar
    4   ${MZ}                # MZ(pole)
    6   ${MT}                # mtop(pole)
Block MINPAR
    1   ${MS}                # Ms
    2   ${Xt}                # Xtt
    3   ${TB}                # TanBeta
Block SPhenoInput       # SPheno specific input 
    1  -1               # error level 
    2   0               # SPA conventions 
    7   0               # Skip 2-loop Higgs corrections 
    8   3               # Method used for two-loop calculation 
    9   1               # Gaugeless limit used at two-loop 
   10   0               # safe-mode used at two-loop 
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
EOF

    rm -f SPheno.spc

    # run the spectrum generator
    slha_output=$($SG SPheno.in SPheno.spc 2>/dev/null)

    if [ -e SPheno.spc ] ; then
        block=$(awk -v block="$output_block" "$print_slha_block_awk" SPheno.spc)
        value=$(echo "$block" | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)
    fi

    [ "x$value" = "x" ] && value="-"

    echo $value
}

run_ss() {
    local SG="$1"
    local MS2=$(echo "scale=5; ${MS}^2" | bc)
    local At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    local Ab=$(echo "scale=10; ${TB} * ${MS}" | bc)
    local Atau=$(echo "scale=10; ${TB} * ${MS}" | bc)
    local MA="$MS"
    local slha_output=
    local block=
    local value=
    local slha_input=
    local output_block=$(echo "${output}" | cut -d'-' -f1)
    local output_entry=$(echo "${output}" | cut -d'-' -f2)

    slha_input=$(
    { cat <<EOF
Block MODSEL                 # Select model
    1    0                   # mSUGRA
   12    ${MS}
Block SOFTSUSY               # SOFTSUSY specific inputs
    1   1.000000000e-05      # tolerance
    2   2.000000000e+00      # up-quark mixing (=1) or down (=2)
    5   1.000000000E+00      # 2-loop running
    3   0.000000000E+00      # printout
    7   2                    # Mh-loops
Block MINPAR                 # Input parameters
    4   1.000000000e+00      # sign(mu)
${sminputs_tmpl}
Block SMINPUTS               # Standard Model inputs
    1   ${AI}                # alpha_em(MZ) SM MSbar
    2   ${GF}                # G_Fermi
    3   ${AS}                # alpha_s(MZ) SM MSbar
    4   ${MZ}                # MZ(pole)
    6   ${MT}                # mtop(pole)
BLOCK EXTPAR
         0     ${MS}   # Q
         1     ${MS}   # M1
         2     ${MS}   # M2
         3     ${M3}   # M3
        11     ${At}   # At
        12     ${Ab}   # Ab
        13     ${Atau} # Atau
        23     ${MS}   # MUE
        26     ${MA}   # MA0
        25     ${TB}   # TB
        31     ${MS}   # MSL(1)
        32     ${MS}   # MSL(2)
        33     ${MS}   # MSL(3)
        34     ${MS}   # MSE(1)
        35     ${MS}   # MSE(2)
        36     ${MS}   # MSE(3)
        41     ${MS}   # MSQ(1)
        42     ${MS}   # MSQ(2)
        43     ${MS}   # MSQ(3)
        44     ${MS}   # MSU(1)
        45     ${MS}   # MSU(2)
        46     ${MS}   # MSU(3)
        47     ${MS}   # MSD(1)
        48     ${MS}   # MSD(2)
        49     ${MS}   # MSD(3)
EOF
    })

    # run the SOFTSUSY spectrum generator
    slha_output=$(echo "$slha_input" | $SG leshouches 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    [ "x$dump_ss_slha_input_file" != "x" ] && \
        echo "$slha_input" > "$dump_ss_slha_input_file"

    [ "x$dump_ss_slha_output_file" != "x" ] && \
        echo "$slha_output" > "$dump_ss_slha_output_file"

    echo $value
}

run_fh() {
    local fh_dir="$1"
    ./run_FeynHiggs_SLHA.sh "$fh_dir" "${MS}" "${TB}" "${Xt}"
}

run_susyhd() {
    local SHDout=
    local Mh=
    local DMh=

    SHDout=$(math -run "Q=${MS}; TB=${TB}; Xt=${Xt}; Get[\"run_SUSYHD.m\"];" 2>&1 >/dev/null)
    Mh=$(echo "$SHDout" | awk '{ print $1 }')
    DMh=$(echo "$SHDout" | awk '{ print $2 }')

    if [ "x$Mh" = "x0" ] ; then
        echo "-   -"
    else
        echo "$Mh   $DMh"
    fi
}

help() {
    cat <<EOF
Usage: $0 [options]
Options:
  --dump-flexiblesusy-slha-input=   dump FlexibleSUSY SLHA input file
  --dump-flexiblesusy-slha-output=  dump FlexibleSUSY SLHA output file
  --dump-softsusy-slha-input=       dump SOFTSUSY SLHA input file
  --dump-softsusy-slha-output=      dump SOFTSUSY SLHA output file
  --output=      output parameter in the format BLOCK-ENTRY1[:ENTRY2] (default: ${output})
  --parameter=   scanned parameter (default: ${parameter})
  --start=       start value (default: ${start})
  --stop=        end value (default: ${stop})
  --steps=       number of steps (default: ${steps})
  --step_size=   linear or log (default: ${step_size})
  --AI=          alpha_em (default: ${AI})
  --AS=          alpha_s (default: ${AS})
  --BL=          beta-functions loop order (default: ${BL})
  --GF           Fermi constant
  --M3factor=    Gluino mass factor: M3 = M3factor * MS (default: ${M3factor})
  --MS=          M_SUSY (default: ${MS})
  --MT=          Top quark pole mass (default: ${MT})
  --MZ           Z pole mass
  --TB=          tan(beta) (default: ${TB})
  --Xt=          Xt (default: ${Xt})
  --help=|-h     print this help message
EOF
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --dump-flexiblesusy-slha-input=*)  dump_fs_slha_input_file=$optarg ;;
            --dump-flexiblesusy-slha-output=*) dump_fs_slha_output_file=$optarg ;;
            --dump-softsusy-slha-input=*)      dump_ss_slha_input_file=$optarg ;;
            --dump-softsusy-slha-output=*)     dump_ss_slha_output_file=$optarg ;;
            --output=*)              output=$optarg ;;
            --parameter=*)           parameter=$optarg ;;
            --start=*)               start=$optarg ;;
            --stop=*)                stop=$optarg ;;
            --steps=*)               steps=$optarg ;;
            --step-size=*)           step_size=$optarg ;;
            --AI=*)                  AI=$optarg ;;
            --AS=*)                  AS=$optarg ;;
            --BL=*)                  BL=$optarg ;;
            --GF=*)                  GF=$optarg ;;
            --M3factor=*)            M3factor=$optarg ;;
            --MS=*)                  MS=$optarg ;;
            --MT=*)                  MT=$optarg ;;
            --MZ=*)                  MZ=$optarg ;;
            --TB=*)                  TB=$optarg ;;
            --Xt=*)                  Xt=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

printf "# MS = ${MS}, TanBeta = ${TB}, Xt = ${Xt}\n"
printf "# %14s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n" "$parameter" "MSSMtower" "EFTtower" "MSSMMuBMu" "HSSUSY" "Softsusy" "MSSMMuBMuSPheno" "FeynHiggs" "DeltaFeynHiggs" "SUSYHD" "DeltaSUSYHD" "SPheno" "SPheno FS-like" "MSSMMuBMuYuatMS" "MSSMMuBMuYuatMSSPheno" "MSSMtower(yt(MS)^0L)"

for i in `seq 0 $steps`; do
    # calculate current value for the scanned variable
    case "$step_size" in
        linear)
            value=$(cat <<EOF | bc
scale=10
$start + ($stop - $start)*${i} / $steps
EOF
                 ) ;;
        log)
            value=$(cat <<EOF | bc -l
scale=10
e(l($start) + (l($stop) - l($start))*${i} / $steps)
EOF
                 ) ;;
        *) echo "Error: unknown step size: $step_size"
           exit 1 ;;
    esac

    eval "${parameter}=${value}"

    M3=$(cat <<EOF | bc
scale=10
$MS * $M3factor
EOF
         )

    # run the spectrum generators
    WRITE_EFT=0
    MTmethod=0
    MF_TL_MATCHING=0
    MhMSSMtower=$(run_sg "models/MSSMtower/run_MSSMtower.x")

    WRITE_EFT=0
    MTmethod=0
    MF_TL_MATCHING=1
    MhMSSMtowerTLMF=$(run_sg "models/MSSMtower/run_MSSMtower.x")

    WRITE_EFT=1
    MTmethod=0
    MF_TL_MATCHING=0
    MhEFTtower=$(run_sg "models/MSSMtower/run_MSSMtower.x")

    WRITE_EFT=0
    MTmethod=0
    MF_TL_MATCHING=0
    MhMSSMMuBMu=$(run_sg "models/MSSMMuBMu/run_MSSMMuBMu.x")

    WRITE_EFT=0
    MTmethod=1
    MF_TL_MATCHING=0
    MhMSSMMuBMuSPheno=$(run_sg "models/MSSMMuBMu/run_MSSMMuBMu.x")

    WRITE_EFT=0
    MTmethod=0
    MF_TL_MATCHING=0
    MhMSSMMuBMuYuatMS=$(run_sg "models/MSSMMuBMuYuatMS/run_MSSMMuBMuYuatMS.x")

    WRITE_EFT=0
    MTmethod=1
    MF_TL_MATCHING=0
    MhMSSMMuBMuYuatMSSPheno=$(run_sg "models/MSSMMuBMuYuatMS/run_MSSMMuBMuYuatMS.x")

    WRITE_EFT=0
    MTmethod=0
    MF_TL_MATCHING=0
    MhHSSUSY=$(run_sg "models/HSSUSY/run_HSSUSY.x")

    MhSoftsusy=$(run_ss "${HOME}/packages/softsusy-3.6.2/softpoint.x")

    FHout=$(run_fh "${HOME}/packages/FeynHiggs-2.11.3/build")
    MhFH=$(echo "$FHout" | awk '{ print $1 }')
    DeltaMhFH=$(echo "$FHout" | awk '{ print $2 }')

    SUSYHDout=$(run_susyhd)
    MhSUSYHD=$(echo "$SUSYHDout" | awk '{ print $1 }')
    DeltaMhSUSYHD=$(echo "$SUSYHDout" | awk '{ print $2 }')

    MhSPheno=$(run_spheno "./SPhenoMSSM")
    MhSPhenoHacked=$(run_spheno "./SPhenoMSSM_FlexibleSUSY_like")

    printf "%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n" "$value" "$MhMSSMtower" "$MhEFTtower" "$MhMSSMMuBMu" "$MhHSSUSY" "$MhSoftsusy" "$MhMSSMMuBMuSPheno" "$MhFH" "$DeltaMhFH" "$MhSUSYHD" "$DeltaMhSUSYHD" "$MhSPheno" "$MhSPhenoHacked" "$MhMSSMMuBMuYuatMS" "$MhMSSMMuBMuYuatMSSPheno" "$MhMSSMtowerTLMF"

done
