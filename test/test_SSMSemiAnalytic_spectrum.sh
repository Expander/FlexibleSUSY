#!/bin/sh

# This script generates an SSM spectrum using both the semi-analytic
# and two-scale solvers, and compares that the resulting spectrums are
# the same.

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
print_block="$UTILSDIR/print_slha_block.awk"
print_block_entry="$UTILSDIR/print_slha_block_entry.awk"
remove_block="$UTILSDIR/remove_slha_block.awk"

semi_analytic_slha="$BASEDIR/../model_files/SSMSemiAnalytic/LesHouches.in.SSMSemiAnalytic"
semi_analytic_input="$BASEDIR/test_SSMSemiAnalytic_spectrum.in.spc"
semi_analytic_output="$BASEDIR/test_SSMSemiAnalytic_spectrum.out.spc"
two_scale_slha="$BASEDIR/../model_files/SSM/LesHouches.in.SSM"
two_scale_input="$BASEDIR/test_SSMSemiAnalytic_spectrum_SSM.in.spc"
two_scale_output="$BASEDIR/test_SSMSemiAnalytic_spectrum_SSM.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
semi_analytic_exe="$BASEDIR/../models/SSMSemiAnalytic/run_SSMSemiAnalytic.x"
two_scale_exe="$BASEDIR/../models/SSM/run_SSM.x"

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

if test ! -x "$semi_analytic_exe"; then
    echo "Error: SSMSemiAnalytic spectrum generator not found: $semi_analytic_exe"
    exit 1
fi

if test ! -x "$two_scale_exe"; then
    echo "Error: SSM spectrum generator not found: $two_scale_exe"
    exit 1
fi

remove_mixing_and_input_blocks() {
    $awk_cmd -f "$remove_block" -v block=UULMIX \
        | $awk_cmd -f "$remove_block" -v block=UDLMIX \
        | $awk_cmd -f "$remove_block" -v block=UURMIX \
        | $awk_cmd -f "$remove_block" -v block=UDRMIX \
        | $awk_cmd -f "$remove_block" -v block=UELMIX \
        | $awk_cmd -f "$remove_block" -v block=UERMIX \
        | $awk_cmd -f "$remove_block" -v block=SCALARMIX \
        | $awk_cmd -f "$remove_block" -v block=SMINPUTS \
        | $awk_cmd -f "$remove_block" -v block=MINPAR \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
        | $awk_cmd -f "$remove_block" -v block=SPINFO \
        | $awk_cmd -f "$remove_block" -v block=MODSEL \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSY
}

# generate point using each solver
cp "$semi_analytic_slha" "$semi_analytic_input"
echo -n "running semi-analytic solver ..."
$semi_analytic_exe --slha-input-file="$semi_analytic_input" \
                   --slha-output-file="$semi_analytic_output"
echo "done"
echo "SSMSemiAnalytic SLHA input file: $semi_analytic_input"
echo "SSMSemiAnalytic SLHA output file: $semi_analytic_output"

if test ! -r "${semi_analytic_output}"; then
    echo "Error: generated SSMSemiAnalytic spectrum not found: $semi_analytic_output"
    exit 1
fi

cp "$two_scale_slha" "$two_scale_input"
echo -n "running two-scale solver ..."
$two_scale_exe --slha-input-file="$two_scale_input" \
               --slha-output-file="$two_scale_output"
echo "done"
echo "SSM SLHA input file: $two_scale_input"
echo "SSM SLHA output file: $two_scale_output"

if test ! -r "${two_scale_output}"; then
    echo "Error: generated SSM spectrum not found: $two_scale_output"
    exit 1
fi

# remove comments from the output files
$sed_cmd -i~ -e '/^ *#/d' "$semi_analytic_output"
$sed_cmd -i~ -e '/^ *#/d' "$two_scale_output"

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities, and input blocks
cp "$semi_analytic_output" "$semi_analytic_output~"
cp "$two_scale_output" "$two_scale_output~"

remove_mixing_and_input_blocks < "$semi_analytic_output~" > "$semi_analytic_output"
remove_mixing_and_input_blocks < "$two_scale_output~" > "$two_scale_output"

echo ""
echo "comparing files $semi_analytic_output and $two_scale_output ..."
echo "required relative agreement: $rel_error"

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 "$semi_analytic_output" "$two_scale_output"`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $semi_analytic_output and  $two_scale_output larger than $rel_error"
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
