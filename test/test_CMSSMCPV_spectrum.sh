#!/bin/bash

# spectrums are the same.

BASEDIR=$(dirname $0)

cmssmcpv_input="$BASEDIR/../model_files/CMSSMCPV/LesHouches.in.CMSSMCPV.wrong_higgs_state"
cmssmcpv_output="$BASEDIR/CMSSMCPV.out.spc"
cmssm_output="$BASEDIR/CMSSM.out.spc"
rel_error="0.0000001"

cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
cmssmcpv_exe="$BASEDIR/../models/CMSSMCPV/run_CMSSMCPV.x"

if test ! -x "$cmssmcpv_exe"; then
    echo "Error: CMSSMCPV spectrum generator not found: $cmssmcpv_exe"
    exit 1
fi

if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi

echo -n "running CMSSMCPV point ... "
$cmssmcpv_exe --slha-input-file=$cmssmcpv_input --slha-output-file=$cmssmcpv_output > /dev/null 2>&1
echo "done"
echo "CMSSMCPV SLHA input file:  $cmssmcpv_input"
echo "CMSSMCPV SLHA output file: $cmssmcpv_output"

echo -n "running CMSSMCPV point ... "
$cmssm_exe --slha-input-file=$cmssmcpv_input --slha-output-file=$cmssm_output > /dev/null 2>&1
echo "done"
echo "CMSSM SLHA input file:  $cmssmcpv_input"
echo "CMSSM SLHA output file: $cmssm_output"

mh_cmssmcpv="`grep hh\(2\) $cmssmcpv_output | awk '{ print $2 }'`"
mh_cmssm="`grep hh\(1\) $cmssm_output | awk '{ print $2 }'`"

# convert scientific notation to bc friendly notation
mh_cmssmcpv="`echo ${mh_cmssmcpv} | sed -e 's/[eE]+*/\\*10\\^/'`"
mh_cmssm="`echo ${mh_cmssm} | sed -e 's/[eE]+*/\\*10\\^/'`"

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=10
abs($mh_cmssm - $mh_cmssmcpv) / ($mh_cmssm) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $mh_cmssmcpv and $mh_cmssm is larger than $rel_error"
    echo "Test status: FAIL"
    exit 1
else
    echo "Test status: OK"
fi

exit 0
