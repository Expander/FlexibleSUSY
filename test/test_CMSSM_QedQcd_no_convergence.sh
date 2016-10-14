#!/bin/sh

BASEDIR=$(dirname $0)
cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
cmssm_in="$BASEDIR/test_CMSSM_QedQcd_no_convergence.in.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi

slha_out=$($cmssm_exe --slha-input-file="$cmssm_in" --slha-output-file=-)

spinfo_4=$(echo "$slha_out" | \
                  awk -f "$print_block" -v block=SPINFO | \
                  awk '{ if ($1 == 4) { $1 = ""; print $0 } }')

error=1

case "$spinfo_4" in
    *no\ convergence*) error=0 ;;
esac

if [ $error -eq 0 ] ; then
    echo "Test result: OK"
else
    echo "Error: point has converged"
    echo ""
    echo "Test result: FAIL"
fi

exit ${error}
