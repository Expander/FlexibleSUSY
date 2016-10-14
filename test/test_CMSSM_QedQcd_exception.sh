#!/bin/sh

BASEDIR=$(dirname $0)
cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
cmssm_in="$BASEDIR/test_CMSSM_QedQcd_exception.in.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi

slha_out=$($cmssm_exe --slha-input-file="$cmssm_in" --slha-output-file=-)

spinfo_4=$(echo "$slha_out" | \
                  awk -f "$print_block" -v block=SPINFO | \
                  awk '{ if ($1 == 4) { $1 = ""; print $0 } }')

error=0

[ -z "$spinfo_4" ] && error=1

if [ $error -eq 0 ] ; then
    echo "Test result: OK"
else
    echo "Error: no SPINFO[4] output entry"
    echo ""
    echo "Test result: FAIL"
fi

exit ${error}
