#!/bin/sh

./models/MSSMNoFV/run_MSSMNoFV.x \
    --slha-input-file=tc.in \
    --slha-output-file=tc.out \
    --rgflow-output-file=tc_rgflow_MSSMNoFV.dat > /dev/null

gnuplot -e "filename='tc_rgflow_MSSMNoFV.dat'; g1column=30; g2column=31" tc_plot_rgflow.gnuplot

./models/phdE6SSM/run_phdE6SSM.x \
    --slha-input-file=tc.in \
    --slha-output-file=tc.out \
    --rgflow-output-file=tc_rgflow_phdE6SSM.dat > /dev/null

gnuplot -e "filename='tc_rgflow_phdE6SSM.dat'; g1column=44; g2column=45" tc_plot_rgflow.gnuplot
