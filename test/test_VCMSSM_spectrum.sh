#!/bin/sh

# This script generates a VCMSSM spectrum and writes it to a SLHA file.
# The EWSB output parameters from this point are then used as inputs for
# the CMSSM.  After running the CMSSM we compare that the VCMSSM and CMSSM
# spectrums are the same.

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
print_block="$UTILSDIR/print_slha_block.awk"
print_block_entry="$UTILSDIR/print_slha_block_entry.awk"
remove_block="$UTILSDIR/remove_slha_block.awk"

vcmssm_slha="$BASEDIR/../model_files/VCMSSM/LesHouches.in.VCMSSM"
vcmssm_input="$BASEDIR/test_VCMSSM_spectrum.in.spc"
vcmssm_output="$BASEDIR/test_VCMSSM_spectrum.out.spc"
cmssm_input="$BASEDIR/test_VCMSSM_spectrum_CMSSM.in.spc"
cmssm_output="$BASEDIR/test_VCMSSM_spectrum_CMSSM.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
vcmssm_exe="$BASEDIR/../models/VCMSSM/run_VCMSSM.x"
cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"

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

if test ! -x "$vcmssm_exe"; then
    echo "Error: VCMSSM spectrum generator not found: $vcmssm_exe"
    exit 1
fi

if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi

remove_mixing_and_input_blocks() {
    $awk_cmd -f "$remove_block" -v block=UMIX \
        | $awk_cmd -f "$remove_block" -v block=VMIX \
        | $awk_cmd -f "$remove_block" -v block=PSEUDOSCALARMIX \
        | $awk_cmd -f "$remove_block" -v block=DSQMIX \
        | $awk_cmd -f "$remove_block" -v block=SELMIX \
        | $awk_cmd -f "$remove_block" -v block=SCALARMIX \
        | $awk_cmd -f "$remove_block" -v block=NMIX \
        | $awk_cmd -f "$remove_block" -v block=CHARGEMIX \
        | $awk_cmd -f "$remove_block" -v block=USQMIX \
        | $awk_cmd -f "$remove_block" -v block=SNUMIX \
        | $awk_cmd -f "$remove_block" -v block=SMINPUTS \
        | $awk_cmd -f "$remove_block" -v block=MINPAR \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
        | $awk_cmd -f "$remove_block" -v block=SPINFO
}

set_parameter_output_scale() {
    if test $# -lt 2; then
        echo "set_parameter_output_scale: Too few arguments"
        return 1
    fi
    $awk_cmd -f "$remove_block" -v block=MODSEL "$1" > "$2"
    _output_scale=$(awk -f "$print_block" -v block=SMINPUTS "$2" \
                        | awk -f "$print_block_entry" -v entries=4)
    echo "Block MODSEL
   12   $_output_scale
" >> $2
}

match_input_blocks() {
    if test $# -lt 3; then
        echo "replace_minpar_block: Too few arguments"
        return 1
    fi
    _input_file="$1"
    _output_file="$2"
    _tb_val="$3"

    $awk_cmd -f "$remove_block" -v block=MINPAR "$_input_file" \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR > "$_output_file"
    $awk_cmd -f "$print_block" -v block=MINPAR "$_input_file" >> "$_output_file"
    echo "    3   $_tb_val" >> "$_output_file"
}

# overwrite parameter output scale
set_parameter_output_scale $vcmssm_slha $vcmssm_input

# generate VCMSSM point
echo -n "running VCMSSM point ..."
$vcmssm_exe --slha-input-file=$vcmssm_input --slha-output-file=$vcmssm_output
echo "done"
echo "VCMSSM SLHA input file:  $vcmssm_input"
echo "VCMSSM SLHA output file: $vcmssm_output"

if test ! -r "$vcmssm_output"; then
    echo "Error: generated VCMSSM spectrum not found: $vcmssm_output"
    exit 1
fi

# remove comments from VCMSSM output spectrum
$sed_cmd -i~ -e '/^ *#/d' $vcmssm_output

# construct matching CMSSM input point
tb_sol=$(awk -f "$print_block" -v block=HMIX "$vcmssm_output" \
             | awk -f "$print_block_entry" -v entries=2)
match_input_blocks $vcmssm_input $cmssm_input $tb_sol

# generate CMSSM point
echo -n "running CMSSM point ..."
$cmssm_exe --slha-input-file=$cmssm_input --slha-output-file=$cmssm_output
echo "done"
echo "CMSSM SLHA input file:  $cmssm_input"
echo "CMSSM SLHA output file: $cmssm_output"

# remove comments from CMSSM output spectrum
$sed_cmd -i~ -e '/^ *#/d' $cmssm_output

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities, and input blocks since these do not match
cp $vcmssm_output $vcmssm_output~
cp $cmssm_output $cmssm_output~

remove_mixing_and_input_blocks < $vcmssm_output~ > $vcmssm_output
remove_mixing_and_input_blocks < $cmssm_output~ > $cmssm_output

echo ""
echo "comparing files $vcmssm_output and $cmssm_output ..."
echo "required relative agreement: $rel_error"

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 $vcmssm_output $cmssm_output`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $vcmssm_output and $cmssm_output larger than $rel_error"
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
