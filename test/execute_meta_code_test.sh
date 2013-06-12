#!/bin/sh

# This script executes all Mathematica tests that are given as
# parameters

if test $# -lt 1 ; then
    echo "Error: Too few arguments"
    echo "./`basename $0` <math-cmd> test.m"
    exit 1
fi

math_cmd="$1"
shift

for t in "$@"
do
    file_ext=$(echo $t | awk -F . '{if (NF>1) {print $NF}}')
    if test "x$file_ext" != "xm"; then
	echo "Error: file does not seem to be a Mathematica script: $t"
        echo "  because the file extension (.$file_ext) is not .m"
	continue;
    fi

    log_file=`echo $t | sed -e 's/\.m$/.log/'`
    rm -f $log_file

    echo -n "executing test: $t ... ";
    echo "**************************************************" >> $log_file;
    echo "* executing test: $t " >> $log_file;
    echo "**************************************************" >> $log_file;

    $math_cmd -run "AppendTo[\$Path, \"./meta/\"]; Get[\"$t\"]; Quit[TestSuite\`GetNumberOfFailedTests[]]" >> $log_file 2>> $log_file;

    if [ $? = 0 ]; then
	echo "OK";
    else
	echo "FAILED";
    fi
done
