#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
FSCONFIG="${HOMEDIR}/flexiblesusy-config"

DEFAULT_CMSSM_INPUT="${HOMEDIR}/model_files/CMSSM/LesHouches.in.CMSSM"
DEFAULT_MSSM_INPUT="${HOMEDIR}/model_files/MSSM/LesHouches.in.MSSM"
DEFAULT_SM_INPUT="${HOMEDIR}/model_files/SM/LesHouches.in.SM"

SGS="\
CMSSM,${DEFAULT_CMSSM_INPUT},0
CMSSMConvergenceTester,${DEFAULT_CMSSM_INPUT},0
CMSSMFPIAbsolute,${DEFAULT_CMSSM_INPUT},0
CMSSMFPIRelative,${DEFAULT_CMSSM_INPUT},0
CMSSMFPITadpole,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLBroyden,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLHybrid,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLHybridS,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLNewton,${DEFAULT_CMSSM_INPUT},0
CMSSMMassWInput,${DEFAULT_CMSSM_INPUT},0
CMSSMNoFV,_DEFAULT_,0
CMSSMCKM,_DEFAULT_,0
CMSSMCPV,_DEFAULT_,0
cCMSSM,_DEFAULT_,0
E6SSM,_DEFAULT_,0
InertMSSM,${DEFAULT_CMSSM_INPUT},0
LHInputMSSM,_DEFAULT_,0
lowMSSM,_DEFAULT_,0
lowNMSSM,_DEFAULT_,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP1_1loop.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP1_2loop_alphaS.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP1_2loop_full.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP2_1loop.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP2_2loop_alphaS.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/TP2_2loop_full.in,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.lowNMSSM_goldstone_tachyon,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.lowNMSSM.pseudoscalar,0
minMSSM,${DEFAULT_CMSSM_INPUT},1
MRSSM,_DEFAULT_,0
MSSM,_DEFAULT_,0
MSSMatMGUT,_DEFAULT_,0
MSSMNoFV,_DEFAULT_,0
MSSMNoFVatMGUT,_DEFAULT_,0
cMSSM,_DEFAULT_,0
cMSSM,${DEFAULT_MSSM_INPUT},0
munuSSM,_DEFAULT_,0
NMSSM,_DEFAULT_,0
NMSSMCPV,_DEFAULT_,0
NoYukawaMSSM,${DEFAULT_CMSSM_INPUT},1
NSM,_DEFAULT_,0
NUHMSSM,_DEFAULT_,0
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP1,1
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP2,1
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP3,0
rootMSSM,${DEFAULT_CMSSM_INPUT},1
SM,_DEFAULT_,0
cSM,${DEFAULT_SM_INPUT},0
SMHighPrecision,${DEFAULT_SM_INPUT},0
SMSSM,_DEFAULT_,0
THDMII,_DEFAULT_,0
TMSSM,_DEFAULT_,0
UMSSM,_DEFAULT_,0
YukawaCMSSM,${DEFAULT_CMSSM_INPUT},0
"

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
    echo "> input file: ${input}"
    echo "> command: ${cmd}"
    eval "${cmd}"

    exit_code="$?"
    echo "> exit code: ${exit_code}"
    echo "> expected result: ${expected_result}"

    if test ${expected_result} -eq ${exit_code} ; then
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
