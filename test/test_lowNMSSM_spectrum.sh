#!/bin/bash

# spectrums are the same.

BASEDIR=$(dirname $0)

lowmssm_input="$BASEDIR/../model_files/lowNMSSM/LesHouches.in.lowNMSSM.pseudoscalar"
lowmssm_output="$BASEDIR/lowNMSSM.out.spc"
softsusy_output="$BASEDIR/SoftsusyNMSSM.out.spc"
rel_error="0.00001"

softsusy_exe="$BASEDIR/../models/SoftsusyNMSSM/run_softpoint.x"
lowmssm_exe="$BASEDIR/../models/lowNMSSM/run_lowNMSSM.x"

if test ! -x "$lowmssm_exe"; then
    echo "Error: lowNMSSM spectrum generator not found: $lowmssm_exe"
    exit 1
fi

if test ! -x "$softsusy_exe"; then
    echo "Error: SoftsusyNMSSM spectrum generator not found: $softsusy_exe"
    exit 1
fi

echo -n "running lowNMSSM point ... "
$lowmssm_exe --slha-input-file=$lowmssm_input --slha-output-file=$lowmssm_output > /dev/null 2>&1
echo "done"
echo "lowNMSSM SLHA input file:  $lowmssm_input"
echo "lowNMSSM SLHA output file: $lowmssm_output"

echo -n "running lowNMSSM point ... "
$softsusy_exe leshouches < $lowmssm_input > $softsusy_output 2> /dev/null
echo "done"
echo "SoftsusyNMSSM SLHA input file:  $lowmssm_input"
echo "SoftsusyNMSSM SLHA output file: $softsusy_output"

mh_lownmssm="`grep hh\(1\) $lowmssm_output | awk '{ print $2 }'`"
mh_softsusy="`grep h0\(1\) $softsusy_output | awk '{ print $2 }'`"

# convert scientific notation to bc friendly notation
mh_lownmssm="`echo ${mh_lownmssm} | sed -e 's/[eE]+*/\\*10\\^/'`"
mh_softsusy="`echo ${mh_softsusy} | sed -e 's/[eE]+*/\\*10\\^/'`"

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=10
abs($mh_softsusy - $mh_lownmssm) / ($mh_softsusy) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $mh_lownmssm and $mh_softsusy is larger than $rel_error"
    echo "Test status: FAIL"
    exit 1
else
    echo "Test status: OK"
fi

exit 0
