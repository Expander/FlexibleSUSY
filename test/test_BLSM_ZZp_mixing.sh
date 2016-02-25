#!/bin/sh

# Check that Z pole mass is close to running Z mass

BASEDIR=$(dirname $0)
CONFIGDIR="${BASEDIR}/../config"
print_block="${CONFIGDIR}/print_slha_block.awk"

blsm_input="$BASEDIR/test_BLSM_ZZp_mixing.in.spc"
blsm_output="$BASEDIR/test_BLSM_ZZp_mixing.out.spc"
blsm_exe="$BASEDIR/../models/BLSM/run_BLSM.x"

if test ! -x "$blsm_exe"; then
    echo "Error: BLSM spectrum generator not found: $blsm_exe"
    exit 1
fi

echo -n "running BLSM point ... "
$blsm_exe --slha-input-file=$blsm_input --slha-output-file=$blsm_output > /dev/null 2>&1
echo "done"
echo "BLSM SLHA input file:  $blsm_input"
echo "BLSM SLHA output file: $blsm_output"

MZin=$(awk -f "$print_block" -v block=SMINPUTS "$blsm_output" | awk '{ if ($1 == 4 ) print $2 }')
MZ=$(  awk -f "$print_block" -v block=MASS "$blsm_output"     | awk '{ if ($1 == 23) print $2 }')
MZp=$( awk -f "$print_block" -v block=MASS "$blsm_output"     | awk '{ if ($1 == 31) print $2 }')

echo ""
echo "input MZ       = $MZin"
echo "calculated MZ  = $MZ"
echo "calculated MZp = $MZp"
echo ""

# convert scientific notation to bc friendly notation
MZin=$(echo "${MZin}" | sed -e 's/[eE]+*/\*10\^/')
MZ=$(  echo "${MZ}"   | sed -e 's/[eE]+*/\*10\^/')
MZp=$( echo "${MZp}"  | sed -e 's/[eE]+*/\*10\^/')

diff=$(cat <<EOF | bc ${BASEDIR}/abs.bc
scale=10
abs($MZ - $MZin) < abs($MZp - $MZin)
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: mass of Z' (${MZp})"
    echo "   is closer to input value ${MZin}"
    echo "   than calculated pole mass of Z ($MZ)"
    echo "Test status: FAIL"
    exit 1
else
    echo "Test status: OK"
fi

exit 0
