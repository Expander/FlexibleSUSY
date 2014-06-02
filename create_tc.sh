#!/bin/sh

./higgs-study/scan_MSSM_tc.x \
    --tc-start=-0.1 \
    --tc-stop=0.1 \
    --tc-npoints=80 \
    --tanb-start=0 \
    --tanb-stop=50 \
    --tanb-npoints=10 > higgs-study/data/tc/scan_MSSMNoFV_tc.dat &

./higgs-study/scan_E6SSM_tc.x \
    --tc-start=-0.1 \
    --tc-stop=0.1 \
    --tc-npoints=80 \
    --tanb-start=0 \
    --tanb-stop=50 \
    --tanb-npoints=10 > higgs-study/data/tc/scan_phdE6SSM_tc.dat &

wait

gnuplot -e "filename='higgs-study/data/tc/scan_MSSMNoFV_tc.dat'" higgs-study/plot_tc.gnuplot

gnuplot -e "filename='higgs-study/data/tc/scan_phdE6SSM_tc.dat'" higgs-study/plot_tc.gnuplot
