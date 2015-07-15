set terminal pdf enhanced
input = '/home/bachi/Programme/FlexibleSUSY/markus_tests/sinThetaW_variation.dat'

set output '/home/bachi/Programme/FlexibleSUSY/markus_tests/sinThetaW_dvbsusy.pdf'
set xlabel 'M_0 / GeV (= M_{12} / GeV)'
set logscale x
set ylabel '1 - sin{/Symbol q}_W({/Symbol d}_{VB}^{SUSY} off) / sin{/Symbol q}_W({/Symbol d}_{VB}^{SUSY} on)'
set key off
plot input using 1:2 with points

set output '/home/bachi/Programme/FlexibleSUSY/markus_tests/sinThetaW_2loopsm.pdf'
set xlabel 'M_0 / GeV (= M_{12} / GeV)'
set logscale x
set ylabel '1 - sin{/Symbol q}_W(SM 2loop off) / sin{/Symbol q}_W(SM 2loop on)'
set key off
plot input using 1:3 with points

# set title 'Schau dir das mal an!!'
# set xrange [-0.5:6.5]
# set yrange [0:5]
# set label 'Hallo' at 1,4