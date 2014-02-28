#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

examples_dir=$(readlink -f "${BASEDIR}/../examples")

rm -rf ${BASEDIR}/standalone/
cp -r ${examples_dir}/standalone/ ${BASEDIR}
(cd ${BASEDIR}/standalone && make)
exit_code="$?"
rm -rf ${BASEDIR}/standalone/

exit ${exit_code}
