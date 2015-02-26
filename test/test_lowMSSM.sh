#!/bin/sh

# This script generates a CMSSM spectrum and writes it to a SLHA file.
# This SLHA file is then used as input point for the low-scale MSSM
# (lowMSSM).  Afterwards we compare that the CMSSM and lowMSSM
# spectrums are the same.

BASEDIR=$(dirname $0)
CONFIGDIR=${BASEDIR}/../config

mssm_input="$BASEDIR/../model_files/CMSSM/LesHouches.in.CMSSM"
mssm_output="$BASEDIR/CMSSM.out.spc"
lowmssm_input="$BASEDIR/lowMSSM.in.spc"
lowmssm_output="$BASEDIR/lowMSSM.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
mssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
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
    echo "Error: CMSSM spectrum generator not found: $mssm_exe"
    exit 1
fi
if test ! -x "$lowmssm_exe"; then
    echo "Error: lowMSSM spectrum generator not found: $lowmssm_exe"
    exit 1
fi

remove_input_blocks() {
    $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=SMINPUTS \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MINPAR \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=EXTPAR \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSOFTIN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=HMIXIN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSQ2IN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSL2IN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSU2IN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSD2IN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=MSE2IN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=TeIN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=TuIN \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=TdIN
}

remove_mixing_matrix_blocks() {
    $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=UMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=VMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=PSEUDOSCALARMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=DSQMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=SELMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=SCALARMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=NMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=CHARGEMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=USQMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=SNUMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=FlexibleSUSYOutput -v entry=0
}

# generate CMSSM point
echo -n "running CMSSM point ... "
$mssm_exe --slha-input-file=$mssm_input --slha-output-file=$mssm_output
echo "done"
echo "CMSSM SLHA input file:  $mssm_input"
echo "CMSSM SLHA output file: $mssm_output"

if test ! -r "$mssm_output"; then
    echo "Error: generated CMSSM spectrum not found: $mssm_output"
    exit 1
fi

# remove comments from CMSSM output spectrum
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
echo "lowMSSM SLHA input file:  $lowmssm_input"
echo "lowMSSM SLHA output file: $lowmssm_output"

# remove comments and input blocks from lowMSSM output spectrum
cp $lowmssm_output $lowmssm_output~
$sed_cmd -e '/^ *#/d' < $lowmssm_output~ | remove_input_blocks > $lowmssm_output

if test ! -r "$lowmssm_output"; then
    echo "Error: generated lowMSSM spectrum not found: $lowmssm_output"
    exit 1
fi

# remove input blocks from CMSSM spectrum file
cp $mssm_output $mssm_output~
remove_input_blocks < $mssm_output~ > $mssm_output

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities
cp $mssm_output $mssm_output~
cp $lowmssm_output $lowmssm_output~

remove_mixing_matrix_blocks < $mssm_output~ > $mssm_output
remove_mixing_matrix_blocks < $lowmssm_output~ > $lowmssm_output

echo ""
echo "comparing files $mssm_output and $lowmssm_output ..."
echo "required relative agreement: $rel_error"

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 $mssm_output $lowmssm_output`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $mssm_output and $lowmssm_output larger that $rel_error"
    echo "$diff"
    echo ""
    echo "Test result: FAIL"
    exit_code=1
else
    echo "$diff"
    echo ""
    echo "Test result: OK"
fi

exit $exit_code
