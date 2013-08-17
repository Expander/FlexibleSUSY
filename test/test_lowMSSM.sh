#!/bin/sh

# This script generates a CMSSM spectrum and writes it to a SLHA file.
# This SLHA file is then used as input point for the low-scale MSSM
# (lowMSSM).  Afterwards we compare that the CMSSM and lowMSSM
# spectrums are the same.

BASEDIR=$(dirname $0)

mssm_input="$BASEDIR/../templates/MSSM/LesHouches.in.MSSM"
mssm_output="$BASEDIR/MSSM.out.spc"
lowmssm_input="$BASEDIR/lowMSSM.in.spc"
lowmssm_output="$BASEDIR/lowMSSM.out.spc"
rel_error="1.3e-4"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
mssm_exe="$BASEDIR/../models/MSSM/run_MSSM.x"
lowmssm_exe="$BASEDIR/../models/lowMSSM/run_lowMSSM.x"

if [ -z "$sed_cmd" ]; then
    echo "Error: sed command not found"
    exit 1
fi
if [ -z "$awk_cmd" ]; then
    echo "Error: awk command not found"
    exit 1
fi
if [ -z "$numdiff_cmd" ]; then
    echo "Error: numdiff command not found"
    exit 1
fi
if test ! -x "$mssm_exe"; then
    echo "Error: MSSM spectrum generator not found: $mssm_exe"
    exit 1
fi
if test ! -x "$lowmssm_exe"; then
    echo "Error: lowMSSM spectrum generator not found: $lowmssm_exe"
    exit 1
fi

# generate CMSSM point
echo -n "running CMSSM point ... "
$mssm_exe --slha-input-file=$mssm_input --slha-output-file=$mssm_output
echo "done"

if test ! -r "$mssm_output"; then
    echo "Error: generated MSSM spectrum not found: $mssm_output"
    exit 1
fi

# remove comments from MSSM output spectrum
$sed_cmd -i~ -e '/^ *#/d' $mssm_output
# rename output blocks to be input blocks for lowMSSM
$sed_cmd \
    -e 's/ MSOFT / MSOFTIN /' \
    -e 's/ HMIX / HMIXIN /' \
    -e 's/ MSQ2 / MSQ2IN /' \
    -e 's/ MSL2 / MSL2IN /' \
    -e 's/ MSU2 / MSU2IN /' \
    -e 's/ MSD2 / MSD2IN /' \
    -e 's/ MSE2 / MSE2IN /' \
    -e 's/ Te / TeIN /' \
    -e 's/ Tu / TuIN /' \
    -e 's/ Td / TdIN /' \
    < $mssm_output > $lowmssm_input

# generate lowMSSM point
echo -n "running lowMSSM point ... "
$lowmssm_exe --slha-input-file=$lowmssm_input --slha-output-file=$lowmssm_output
echo "done"

# remove comments from lowMSSM output spectrum
cp $lowmssm_output $lowmssm_output~
$sed_cmd -e '/^ *#/d' < $lowmssm_output~ | $awk_cmd -f $BASEDIR/remove_input_blocks.awk > $lowmssm_output

if test ! -r "$lowmssm_output"; then
    echo "Error: generated lowMSSM spectrum not found: $lowmssm_output"
    exit 1
fi

diff=`$numdiff_cmd --relative-tolerance $rel_error $mssm_output $lowmssm_output`
diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $mssm_output and $lowmssm_output larger that $rel_error"
    echo "$diff"
    exit 1
else
    echo "$diff"
fi
