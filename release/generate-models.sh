#!/bin/sh

# creates models for public release

number_of_jobs=1

#_____________________________________________________________________
help() {
cat <<EOF
Usage: ./`basename $0` [options]
Options:

  --number-of-jobs=    number of parallel makefile jobs
  --help,-h            Print this help message
EOF
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --number-of-jobs=*)      number_of_jobs=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

echo "Using $number_of_jobs parallel jobs"

# models and corresponding SARAH model files
models_array="             \
   CMSSM,MSSM              \
   MSSM,MSSM               \
   MSSMatMGUT,MSSM         \
   MSSMNoFV,MSSMNoFV       \
   MSSMNoFVatMGUT,MSSMNoFV \
   CMSSMNoFV,MSSMNoFV      \
   NUHMSSM,MSSM            \
   lowMSSM,MSSM            \
   MSSMRHN,MSSMRHN         \
   NMSSM,NMSSM             \
   NUTNMSSM,NMSSM          \
   SMSSM,SMSSM             \
   UMSSM,UMSSM             \
   E6SSM,E6SSM             \
   MRSSM,MRSSM             \
   TMSSM,TMSSM             \
"

models="`echo $models_array | sed 's/,[a-zA-Z0-9]*//g'`"

# directory of this script
BASEDIR=$(dirname $0)

# creating models
for m in ${models_array}
do
    _model_sarah="`echo $m | tr ',' ' '`"
    _model="`echo $_model_sarah | cut -d\" \" -f1`"
    _sarah="`echo $_model_sarah | cut -d\" \" -f2`"

    ./createmodel --name=${_model} --sarah-model=${_sarah} --force
done

# running configure
models_comma="`echo $models | tr ' ' ','`"

./configure \
    --disable-compile \
    --with-models=${models_comma}

make showbuild

make -j${number_of_jobs}

# packing models
for m in ${models}
do
    echo "packing ${m} ..."
    make pack-${m}-src
done
