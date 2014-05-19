#!/bin/sh

./models/MSSMNoFV/run_MSSMNoFV.x \
    --slha-input-file=tc.in \
    --slha-output-file=tc.out \
    --rgflow-output-file=tc_rgflow.dat > /dev/null


gnuplot -e "filename='tc_rgflow.dat'" tc_plot_rgflow.gnuplot
