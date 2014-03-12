#!/bin/sh

TIMEFORMAT='%U'
LC_NUMERIC=C
LC_COLLATE=C

# directory of this script
BASEDIR=$(dirname $0)

ce6ssm_specgen="/home/avoigt/tu/research/cE6SSM_SpecGen/softsusy.x"
fs_specgen="${BASEDIR}/models/TestE6SSM/run_TestE6SSM.x"
random_float="${BASEDIR}/random_float.x"
slha_template="${BASEDIR}/ce6ssm_generic.slha2"
slha_file="${BASEDIR}/${slha_template}.point"

out_ce6ssm="out.spc.ce6ssm"
out_fse6ssm="out.spc.fse6ssm"

help() {
cat <<EOF
   --lambda              universal lambda
   --kappa               universal kappa
   --tan_beta            tan(beta)
   --vs                  vs
   --help, -h            help message
EOF
}

# default values
lambda=0.2
lambda12="$lambda"
kappa=0.2
tan_beta=10
vs=1000
muPrime="10000"
BmuPrime="10000"
solution="1"

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --lambda=*)              lambda=$optarg ;;
            --kappa=*)               kappa=$optarg ;;
            --tan_beta=*)            tan_beta=$optarg ;;
            --vs=*)                  vs=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

rm -f $out_ce6ssm $out_fse6ssm

echo "Parameter point: "
printf "%10g " "$lambda"
printf "%10g " "$lambda12"
printf "%10g " "$kappa"
printf "%10g " "$tan_beta"
printf "%10g " "$vs"
printf "%10g " "$muPrime"
printf "%10g " "$BmuPrime"
printf "\n"

echo -n "\nRunning hand-written CE6SSM spectrum generator ... "

# run CE6SSM spectrum generator
error="0"
$ce6ssm_specgen --brief \
    --tan-beta $tan_beta \
    --lambda1 $lambda12 \
    --lambda2 $lambda12 \
    --lambda3 $lambda \
    --kappa1 $kappa \
    --kappa2 $kappa \
    --kappa3 $kappa \
    --vs $vs \
    --mu-prime $muPrime \
    --output-file $out_ce6ssm > /dev/null

error="$?"

if test ! "x${error}" = "x0" ; then
    echo "not ok"
    echo ""
    echo "Error: hand-written CE6SSM spectrum generator could not"
    echo "   calculate spectrum (error code ${error})"
    exit 1
else
    echo "ok"
fi

m0=` awk '{ if ($2 ~ /m0$/ ) print $1 }' $out_ce6ssm`
m12=`awk '{ if ($2 ~ /M12$/) print $1 }' $out_ce6ssm`
a0=` awk '{ if ($2 ~ /A$/  ) print $1 }' $out_ce6ssm`

cp $slha_template $slha_file
echo "Block MINPAR"     >> $slha_file
echo "   1   $m0        # m0"        >> $slha_file
echo "   2   $m12       # m12"       >> $slha_file
echo "   3   $tan_beta  # tan(beta)" >> $slha_file
echo "   5   $a0        # a0"        >> $slha_file
echo "Block EXTPAR"     >> $slha_file
echo "  61   $lambda    # lambda"    >> $slha_file
echo "  62   $kappa     # kappa"     >> $slha_file
echo "  63   $muPrime   # mu'"       >> $slha_file
echo "  64   $BmuPrime  # Bmu'"      >> $slha_file
echo "  65   $vs        # vs"        >> $slha_file
echo "  66   $lambda12  # lambda12"  >> $slha_file

$fs_specgen \
    --slha-input-file=$slha_file \
    --slha-output-file=$out_fse6ssm > /dev/null

error="$?"

if test ! "x${error}" = "x0" ; then
    echo "Error: FlexibleSUSY's CE6SSM spectrum generator could not"
    echo "   calculate spectrum (error code ${error})"
fi

exit $error
