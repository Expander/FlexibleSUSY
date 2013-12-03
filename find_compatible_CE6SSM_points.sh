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

while [ true ]
do
    # create CE6SSM point
    lambda=`$random_float 0.1 0.5`
    kappa=`$random_float 0.1 0.5`
    tan_beta=`$random_float 5 20`
    vs=`$random_float 1000 10000`
    muPrime="10000"
    BmuPrime="10000"
    solution="1"

    rm -f out.spc

    # run CE6SSM spectrum generator
    error="0"
    $ce6ssm_specgen --brief \
        --tan-beta $tan_beta \
        --lambda1 $lambda \
        --lambda2 $lambda \
        --lambda3 $lambda \
        --kappa1 $kappa \
        --kappa2 $kappa \
        --kappa3 $kappa \
        --vs $vs \
        --mu-prime $muPrime \
        --output-file out.spc > /dev/null

    error="$?"

    if test ! "x${error}" = "x0" ; then
        continue
    fi

    m0=` awk '{ if ($2 ~ /m0$/ ) print $1 }' out.spc`
    m12=`awk '{ if ($2 ~ /M12$/) print $1 }' out.spc`
    a0=` awk '{ if ($2 ~ /A$/  ) print $1 }' out.spc`

    slha_file="${BASEDIR}/${slha_template}.point"

    cp $slha_template $slha_file
    echo "Block MINPAR"     >> $slha_file
    echo "   1   $m0"       >> $slha_file
    echo "   2   $m12"      >> $slha_file
    echo "   3   $tan_beta" >> $slha_file
    echo "   5   $a0"       >> $slha_file
    echo "Block EXTPAR"     >> $slha_file
    echo "  61   $lambda"   >> $slha_file
    echo "  62   $kappa"    >> $slha_file
    echo "  63   $muPrime"  >> $slha_file
    echo "  64   $BmuPrime" >> $slha_file
    echo "  65   $vs"       >> $slha_file

    $fs_specgen \
        --slha-input-file=$slha_file \
        --slha-output-file=out.spc > /dev/null

    error="$?"

    if test ! "x${error}" = "x0" ; then
        continue
    fi

    printf "%10g " "$lambda"
    printf "%10g " "$kappa"
    printf "%10g " "$tan_beta"
    printf "%10g " "$vs"
    printf "%10g " "$muPrime"
    printf "%10g " "$BmuPrime"
    printf "\n"
done
