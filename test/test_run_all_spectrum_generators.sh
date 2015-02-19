#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
FSCONFIG="${HOMEDIR}/flexiblesusy-config"

DEFAULT_CMSSM_INPUT="${HOMEDIR}/model_files/CMSSM/LesHouches.in.CMSSM"

SGS="\
CMSSM,${DEFAULT_CMSSM_INPUT},ok
CMSSMConvergenceTester,${DEFAULT_CMSSM_INPUT},ok
CMSSMFPIAbsolute,${DEFAULT_CMSSM_INPUT},ok
CMSSMFPIRelative,${DEFAULT_CMSSM_INPUT},ok
CMSSMFPITadpole,${DEFAULT_CMSSM_INPUT},ok
CMSSMGSLBroyden,${DEFAULT_CMSSM_INPUT},ok
CMSSMGSLHybrid,${DEFAULT_CMSSM_INPUT},ok
CMSSMGSLHybridS,${DEFAULT_CMSSM_INPUT},ok
CMSSMGSLNewton,${DEFAULT_CMSSM_INPUT},ok
CMSSMMassWInput,${DEFAULT_CMSSM_INPUT},ok
CMSSMNoFV,_DEFAULT_,ok
E6SSM,_DEFAULT_,ok
InertMSSM,${DEFAULT_CMSSM_INPUT},ok
LHInputMSSM,_DEFAULT_,ok
lowMSSM,_DEFAULT_,ok
lowNMSSM,_DEFAULT_,ok
minMSSM,${DEFAULT_CMSSM_INPUT},fail
MRSSM,_DEFAULT_,ok
MSSM,_DEFAULT_,ok
MSSMatMGUT,_DEFAULT_,ok
MSSMNoFV,_DEFAULT_,ok
MSSMNoFVatMGUT,_DEFAULT_,ok
munuSSM,_DEFAULT_,ok
NMSSM,_DEFAULT_,ok
NoYukawaMSSM,${DEFAULT_CMSSM_INPUT},fail
NSM,_DEFAULT_,ok
NUHMNMSSM,_DEFAULT_,ok
NUHMSSM,_DEFAULT_,ok
NUHNMSSM,_DEFAULT_,ok
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP1,ok
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP2,fail
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP3,ok
rootMSSM,${DEFAULT_CMSSM_INPUT},fail
SM,_DEFAULT_,ok
SMSSM,_DEFAULT_,ok
THDMII,_DEFAULT_,ok
TMSSM,_DEFAULT_,ok
UMSSM,_DEFAULT_,ok
YukawaCMSSM,${DEFAULT_CMSSM_INPUT},ok
"

[ $("$FSCONFIG" --with-CMSSM) = yes ] &&
    SGS="$SGS customized-betas"

TMP_FILE="${BASEDIR}/test_spectrum_generator.sh.tmp"

rm -f "$TMP_FILE"

errors=0

for setup in ${SGS}
do
    model="`echo ${setup} | tr ',' ' ' | awk '{ print $1 }'`"
    input="`echo ${setup} | tr ',' ' ' | awk '{ print $2 }'`"
    expected_result="`echo ${setup} | tr ',' ' ' | awk '{ print $3 }'`"

    echo "== $model ===================================="

    if [ $("$FSCONFIG" --with-${model}) = no ] ; then
        echo "> skipping, because the model is not configured"
        continue
    fi

    sg="${HOMEDIR}/models/${model}/run_${model}.x"

    if test ! -x "${sg}"; then
        echo "> Error: spectrum generator not built: $sg"
        errors=1
        continue
    fi

    case "$input" in
        _DEFAULT_) input="${HOMEDIR}/model_files/${model}/LesHouches.in.${model}" ;;
    esac

    if test ! -e "${input}"; then
        echo "> Error: input file not found: $input"
        errors=1
        continue
    fi

    cmd="${sg} --slha-input-file=${input} --slha-output-file=${TMP_FILE} > /dev/null 2>&1"

    echo "> running spectrum generator for $model"
    echo "> cmd: ${cmd}"
    eval "${cmd}"

    exit_code="$?"
    echo "> exit code: ${exit_code}"
    echo "> expected result: ${expected_result}"

    if test \( "${expected_result}" = ok   -a ${exit_code} -eq 0 \) -o \
            \( "${expected_result}" = fail -a ${exit_code} -ne 0 \)
    then
        echo "> spectrum generator: ok"
    else
        echo "> spectrum generator: FAIL"
        errors=1
    fi
done

rm -f "$TMP_FILE"

echo ""

if test ${errors} -eq 0 ; then
    echo "Test result: OK"
else
    echo "Test result: FAIL"
fi

exit ${errors}
