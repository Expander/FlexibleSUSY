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

# 7 : MS
# 8 : vu of MSSMtower
# 9 : vu of EFTtower
# 10: vu of MSSMMuBMu
# 11: vu of HSSUSY
# 12: vu of SOFTSUSY

# 13: MS
# 14: v of MSSMtower
# 15: v of EFTtower
# 16: v of MSSMMuBMu
# 17: v of HSSUSY
# 18: v of SOFTSUSY

plot [:] [:] \
     filename u ($1/1000):(calc_mt($2,$8)) t 'MSSM-tower/MSSM' w lines ls 1, \
     filename u ($1/1000):(calc_mt($3,$15)) t 'MSSM-tower/SM' w lines ls 5, \
     "< awk '{ if ($4 > 0) print }' ".filename u ($1/1000):(calc_mt($4,$10)) t 'MSSM 2L' w lines ls 3, \
     filename u ($1/1000):(calc_mt($5,$17)) t 'HSSUSY 2L' w lines ls 2, \
     filename u ($1/1000):(calc_mt($6,$18)) t 'SOFTSUSY 2L' w lines ls 4
