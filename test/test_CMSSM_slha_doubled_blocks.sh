#!/bin/sh

# This checks that doubled SLHA input blocks do not change the output
# spectrum

BASEDIR=$(dirname $0)

mssm_exe="$BASEDIR/../models/MSSM/run_MSSM.x"
mssm_input="$BASEDIR/test_MSSM_slha_doubled_blocks.spc.in"
mssm_output="$BASEDIR/test_MSSM_slha_doubled_blocks.spc.out"

diff_cmd=`command -v diff`

if [ -z "$diff_cmd" ]; then
    echo "Error: diff command not found"
    exit 1
fi

if test ! -x "$mssm_exe"; then
    echo "Error: MSSM spectrum generator not found: $mssm_exe"
    exit 1
fi

echo -n "running CMSSM point ... "
$mssm_exe --slha-input-file=$mssm_input --slha-output-file=$mssm_output
echo "done"
echo "CMSSM SLHA input file:  $mssm_input"
echo "CMSSM SLHA output file: $mssm_output"

# append redundant SMINPUTS block
mssm_input_appended_sminputs="$BASEDIR/test_MSSM_slha_doubled_blocks_appenden_sminputs.spc.in"
mssm_output_appended_sminputs="$BASEDIR/test_MSSM_slha_doubled_blocks_appenden_sminputs.spc.out"

cp $mssm_input $mssm_input_appended_sminputs
echo "Block SMINPUTS" >> $mssm_input_appended_sminputs

echo -n "running CMSSM point with extra appended SMINPUTS block ... "
$mssm_exe --slha-input-file=$mssm_input_appended_sminputs --slha-output-file=$mssm_output_appended_sminputs
echo "done"
echo "CMSSM SLHA input file:  $mssm_input_appended_sminputs"
echo "CMSSM SLHA output file: $mssm_output_appended_sminputs"

difference=`$diff_cmd $mssm_output $mssm_output_appended_sminputs`

if [ -n "$difference" ]; then
    echo "Error: difference between $mssm_output and $mssm_output_appended_sminputs not empty!"
    echo "$difference"
    exit 1
else
    echo "Difference between $mssm_output and $mssm_output_appended_sminputs is empty"
    echo "OK"
fi
