#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

FSCONFIG="$BASEDIR/../flexiblesusy-config"
VALGRIND="`command -v valgrind`"
MSSMtower_EXE="$BASEDIR/../models/MSSMtower/run_MSSMtower.x"
MSSMtower_INPUT="$BASEDIR/../model_files/MSSMtower/LesHouches.in.MSSMtower"
MSSMtower_OUTPUT="$BASEDIR/LesHouches.out.MSSMtower"
VALGRIND_OUTPUT="${BASEDIR}/test_MSSMtower_profile.valgrind.out"
CALLGRIND_OUTPUT="${BASEDIR}/test_MSSMtower_profile.callgrind.out"

if [ -z "$VALGRIND" ]; then
    echo "Error: valgrind not found"
    echo "Profiling will not be performed."
    exit 0
fi

if test ! -x "$MSSMtower_EXE"; then
    echo "Error: MSSMtower spectrum generator not found: $MSSMtower_EXE"
    exit 1
fi

if test ! -f "$MSSMtower_INPUT"; then
    echo "Error: MSSMtower SLHA input file not found: $MSSMtower_INPUT"
    exit 1
fi


valgrind \
    --tool=callgrind \
    --separate-threads=yes \
    --dump-instr=yes \
    --simulate-cache=yes \
    --collect-jumps=yes \
    --log-file=${VALGRIND_OUTPUT} \
    --callgrind-out-file=${CALLGRIND_OUTPUT} \
    ${MSSMtower_EXE} \
    --slha-input-file=${MSSMtower_INPUT} \
    --slha-output-file=${MSSMtower_OUTPUT}

echo "Valgrind log file: ${VALGRIND_OUTPUT}"
echo "Callgring output file: ${CALLGRIND_OUTPUT}"
echo ""
cat ${VALGRIND_OUTPUT}

exit 0
