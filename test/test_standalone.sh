#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

examples_dir=$(readlink -f "${BASEDIR}/../examples")

standalone_dirs="standalone-model standalone-rge"

for dir in ${standalone_dirs}
do
    echo "> cleaning: rm -rf ${BASEDIR}/${dir}"
    rm -rf ${BASEDIR}/${dir}

    echo "> copying: cp -r ${examples_dir}/${dir}/ ${BASEDIR}"
    cp -r ${examples_dir}/${dir}/ ${BASEDIR}

    echo "> building: (cd ${BASEDIR}/${dir} && make)"
    (cd ${BASEDIR}/${dir} && make)

    exit_code="$?"
    echo "> exit code: ${exit_code}"

    if test ${exit_code} -ne 0; then
        echo "=============="
        echo "Result: FAILED"
        echo "Error: build failed in directory ${dir}"
        echo "=============="
        break
    else
        echo "=========="
        echo "Result: OK"
        echo "=========="
    fi

    echo "> cleaning: rm -rf ${BASEDIR}/${dir}"
    rm -rf ${BASEDIR}/${dir}/
done

exit ${exit_code}
