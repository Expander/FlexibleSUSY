#!/bin/sh

# This test compares the real CMSSM with the complex CMSSM (cCMSSM)
# for real input parameters

BASEDIR=$(dirname $0)
CONFIGDIR=${BASEDIR}/../config

cmssm_input="$BASEDIR/../model_files/CMSSM/LesHouches.in.CMSSM"
cmssm_output="$BASEDIR/test_cCMSSM_CMSSM.out.spc"
ccmssm_input="$cmssm_input"
ccmssm_output="$BASEDIR/test_cCMSSM_cCMSSM.out.spc"
ccmssm_real_output="$BASEDIR/test_cCMSSM_cCMSSM_real_part.out.spc"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
ccmssm_exe="$BASEDIR/../models/cCMSSM/run_cCMSSM.x"

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
if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi
if test ! -x "$ccmssm_exe"; then
    echo "Error: cCMSSM spectrum generator not found: $ccmssm_exe"
    exit 1
fi

remove_imaginary_blocks() {
    $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=ImAu \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=ImAd \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=ImAe \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=ImHMIX \
        | $awk_cmd -f $CONFIGDIR/remove_slha_block.awk -v block=ImMSOFT
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
$cmssm_exe --slha-input-file=$cmssm_input --slha-output-file=$cmssm_output
echo "done"
echo "CMSSM SLHA input file:  $cmssm_input"
echo "CMSSM SLHA output file: $cmssm_output"

# generate cCMSSM point
echo -n "running cCMSSM point ... "
$ccmssm_exe --slha-input-file=$ccmssm_input --slha-output-file=$ccmssm_output
echo "done"
echo "cCMSSM SLHA input file:  $ccmssm_input"
echo "cCMSSM SLHA output file: $ccmssm_output"

remove_imaginary_blocks < "$ccmssm_output" > "$ccmssm_real_output"

# remove mixing matrix blocks and comments
cp "$cmssm_output" "$cmssm_output~"
cp "$ccmssm_real_output" "$ccmssm_real_output~"

remove_mixing_matrix_blocks < "$cmssm_output~"       | $sed_cmd -e '/#.*/d' > "$cmssm_output"
remove_mixing_matrix_blocks < "$ccmssm_real_output~" | $sed_cmd -e '/#.*/d' > "$ccmssm_real_output"

echo ""
echo "comparing files: $cmssm_output and $ccmssm_real_output ..."

diff=`$numdiff_cmd\
 $cmssm_output $ccmssm_real_output`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $cmssm_output and $ccmssm_real_output"
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
