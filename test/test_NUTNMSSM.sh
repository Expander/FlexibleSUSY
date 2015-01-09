#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

exit_code=0

for point in BP1 BP2 BP3
do
    echo "Running $point ..."

    input_file="${BASEDIR}/../model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_${point}"
    output_file="${BASEDIR}/LesHouches.out.NUTNMSSM_${point}"

    ${BASEDIR}/../models/NUTNMSSM/run_NUTNMSSM.x \
        --slha-input-file=${input_file} \
        --slha-output-file=${output_file}

    if test "x$?" = "x0"; then
        echo "=========="
        echo "${point}: OK"
        echo "=========="
        grep hh ${output_file}
    else
        echo "=========="
        echo "${point}: FAIL"
        echo "=========="
        grep Problems ${output_file}
        exit_code=1
    fi

    echo ""

    rm -f ${output_file}
done

exit ${exit_code}
