#!/bin/sh

output="MASS-25"
parameter=MS
Xt=0
TB=5
MS=1000
At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
M3factor=1
M3=
AS="1.184000000e-01"
MT="1.733400000e+02"
MTmethod=0
UseMTmethod="$MTmethod"

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
Block FlexibleSUSY
   17   ${UseMTmethod}       # mt calculation (0 = FlexibleSUSY, 1 = SPheno)
Block SMINPUTS               # Standard Model inputs
    3   ${AS}                # alpha_s(MZ) SM MSbar
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
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk")

    [ "x$value" = "x" ] && value="-"

    [ "x$dump_fs_slha_input_file" != "x" ] && \
        echo "$slha_input" > "$dump_fs_slha_input_file"

    [ "x$dump_fs_slha_output_file" != "x" ] && \
        echo "$slha_output" > "$dump_fs_slha_output_file"

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
    3   ${AS}                # alpha_s(MZ) SM MSbar
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
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk")

    [ "x$value" = "x" ] && value="-"

    [ "x$dump_ss_slha_input_file" != "x" ] && \
        echo "$slha_input" > "$dump_ss_slha_input_file"

    [ "x$dump_ss_slha_output_file" != "x" ] && \
        echo "$slha_output" > "$dump_ss_slha_output_file"

    echo $value
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
  --AS=          alpha_s (default: ${AS})
  --M3factor=    Gluino mass factor: M3 = M3factor * MS (default: ${M3factor})
  --MS=          M_SUSY (default: ${MS})
  --MT=          Top quark pole mass (default: ${MT})
  --MTmethod=    0 = FlexibleSUSY, 1 = SPheno (default: $MTmethod)
                 (Only used in FlexibleSUSY/MSSMMuBMu)
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
            --AS=*)                  AS=$optarg ;;
            --M3factor=*)            M3factor=$optarg ;;
            --MS=*)                  MS=$optarg ;;
            --MT=*)                  MT=$optarg ;;
            --MTmethod=*)            MTmethod=$optarg ;;
            --TB=*)                  TB=$optarg ;;
            --Xt=*)                  Xt=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

printf "# MS = ${MS}, TanBeta = ${TB}, Xt = ${Xt}\n"
printf "# %14s %16s %16s %16s %16s\n" "$parameter" "MSSMtower" "MSSMMuBMu" "HSSUSY" "Softsusy"

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
    UseMTmethod=0
    MhMSSMtower=$(run_sg "models/MSSMtower/run_MSSMtower.x")

    UseMTmethod="$MTmethod"
    MhMSSMMuBMu=$(run_sg "models/MSSMMuBMu/run_MSSMMuBMu.x")

    UseMTmethod=0
    MhHSSUSY=$(run_sg "models/HSSUSY/run_HSSUSY.x")

    MhSoftsusy=$(run_ss "${HOME}/packages/softsusy-3.6.2/softpoint.x")

    printf "%16s %16s %16s %16s %16s\n" "$value" "$MhMSSMtower" "$MhMSSMMuBMu" "$MhHSSUSY" "$MhSoftsusy"

done
