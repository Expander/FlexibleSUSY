set terminal pdfcairo
set output "mt_".filename.".pdf"
set key box top right
set logscale x
set grid

set style line 1 lt 1 dt 1 lw 2 lc rgb "#FF0000"
set style line 2 lt 1 dt 2 lw 2 lc rgb "#0000FF"
set style line 3 lt 1 dt 4 lw 2 lc rgb "#45AD53"
set style line 4 lt 1 dt 3 lw 2 lc rgb "#FFBF00"
set style line 5 lt 1 dt 5 lw 2 lc rgb "#FF00FF"

set xlabel "M_S / TeV"
set ylabel "m_t(M_S) / GeV"

calc_mt(yt,v) = yt * v / sqrt(2)

# 1: MS
# 2: yt of MSSMtower
# 3: yt of EFTtower
# 4: yt of MSSMMuBMu
# 5: yt of HSSUSY
# 6: yt of SOFTSUSY
# 7: yt of MSSMMuBMu/SPheno-like

# 8 : MS
# 9 : vu of MSSMtower
# 10: vu of EFTtower
# 11: vu of MSSMMuBMu
# 12: vu of HSSUSY
# 13: vu of SOFTSUSY
# 14: vu of MSSMMuBMu/SPheno-like

# 15: MS
# 16: v of MSSMtower
# 17: v of EFTtower
# 18: v of MSSMMuBMu
# 19: v of HSSUSY
# 20: v of SOFTSUSY
# 21: v of MSSMMuBMu/SPheno-like

plot [:] [:] \
     filename u ($1/1000):(calc_mt($2,$9)) t 'MSSM-tower/MSSM' w lines ls 1, \
     filename u ($1/1000):(calc_mt($3,$17)) t 'MSSM-tower/SM' w lines ls 5, \
     "< awk '{ if ($4 > 0) print }' ".filename u ($1/1000):(calc_mt($4,$11)) t 'MSSM 2L' w lines ls 3, \
     filename u ($1/1000):(calc_mt($5,$19)) t 'HSSUSY 2L' w lines ls 2, \
     filename u ($1/1000):(calc_mt($6,$20)) t 'SOFTSUSY 2L' w lines ls 4
