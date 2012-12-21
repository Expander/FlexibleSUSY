#!/bin/sh

# This script executes all tests that are given as parameters

for t in "$@"
do
    if test ! -x $t; then
	echo "Error: file does not exist or is not executable: $t"
	continue;
    fi

    log_file=`echo $t | sed -e 's/.x$/.log/'`
    rm -f $log_file

    echo -n "executing test: $t ... ";
    echo "**************************************************" >> $log_file;
    echo "* executing test: $t " >> $log_file;
    echo "**************************************************" >> $log_file;

    ./$t --log_level=test_suite >> $log_file 2>> $log_file;

    if [ $? = 0 ]; then
	echo "OK";
    else
	echo "FAILED";
    fi
done
