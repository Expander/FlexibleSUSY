#!/bin/sh

./models/MSSMNoFV/run_MSSMNoFV.x \
    --slha-input-file=tc.in \
    --slha-output-file=tc.out \
    --rgflow-output-file=tc_rgflow_MSSMNoFV.dat > /dev/null

./models/phdE6SSM/run_phdE6SSM.x \
    --slha-input-file=tc.in \
    --slha-output-file=tc.out \
    --rgflow-output-file=tc_rgflow_phdE6SSM.dat > /dev/null

gnuplot -e "filename1='tc_rgflow_MSSMNoFV.dat'; filename2='tc_rgflow_phdE6SSM.dat'; file1g1column=30; file1g2column=31; file2g1column=44; file2g2column=45" tc_plot_rgflow.gnuplot
