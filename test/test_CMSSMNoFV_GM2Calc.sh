#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
GM2CALCDIR=$(readlink -f "${BASEDIR}/../addons/GM2Calc")
FSCONFIG="$BASEDIR/../flexiblesusy-config"
GM2CALC_EXE=$(readlink -f "${BASEDIR}/../addons/GM2Calc/gm2calc.x")
CMSSMNoFV_EXE=$(readlink -f "${BASEDIR}/../models/CMSSMNoFV/run_CMSSMNoFV.x")
SLHA_IN=$(readlink -f "${BASEDIR}/../models/CMSSMNoFV/LesHouches.in.CMSSMNoFV")
SLHA_OUT=$(readlink -f "${BASEDIR}/test_CMSSMNoFV_GM2Calc.out.spc")
print_block="$BASEDIR/../utils/print_slha_block.awk"
rel_error=0.000000001

[ $("$FSCONFIG" --with-CMSSMNoFV) = yes -a -x ${CMSSMNoFV_EXE} ] || {
    echo "Error: CMSSMNoFV needs to be build!"
    exit 1;
}

[ $("$FSCONFIG" --with-GM2Calc) = yes -a -x ${GM2CALC_EXE} ] || {
    echo "Error: GM2Calc needs to be build!"
    exit 1;
}

{ cat "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    15  1   # calculate observables (a_muon, ...)
EOF
  } | "${CMSSMNoFV_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: ${CMSSMNoFV_EXE} failed!"
    exit 1
}

amu_gm2calc_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 1) print $2 }')

amu_gm2calc=$({ cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     4  0  # verbose output
EOF
  cat "${SLHA_OUT}";
      } | "${GM2CALC_EXE}" --slha-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# convert scientific notation to bc friendly notation
amu_gm2calc_fs=$(echo "${amu_gm2calc_fs}" | sed -e 's/[eE]/\*10\^/')
amu_gm2calc=$(echo "${amu_gm2calc}" | sed -e 's/[eE]/\*10\^/')

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=10
abs((abs($amu_gm2calc_fs) - abs($amu_gm2calc)) / ($amu_gm2calc_fs)) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $amu_gm2calc_fs and $amu_gm2calc is larger than $rel_error"
    echo "Test status: FAIL"
    exit 1
else
    echo "FlexibleSUSY: amu = $amu_gm2calc_fs"
    echo "GM2Calc     : amu = $amu_gm2calc"
    echo "Test status : OK"
fi
