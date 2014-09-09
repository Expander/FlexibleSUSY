set terminal pdf

set title "MSSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MSSM_rgflow.dat'
set output "MSSM_rgflow.pdf"

negsqrt(x) = x < 0 ? -sqrt(-x) : sqrt(x)

# parameter  column
# Mu         29
# BMu        62
# vu         34
# vd         35

set style line 1 lt 1 lw 2 lc rgb "red"
set style line 2 lt 2 lw 2 lc rgb "black"
set style line 3 lt 3 lw 2 lc rgb "blue"
set style line 4 lt 4 lw 2 lc rgb "cyan"
set style line 5 lt 5 lw 2 lc rgb "yellow"

plot \
     for [i=2:2]     filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=6:6]     filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=10:10]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=11:11]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=15:15]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=19:19]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=20:20]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=24:24]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=28:28]   filename using 1:(column(i))       title columnhead(i) w lines ls 3, \
     for [i=30:32]   filename using 1:(column(i))       title columnhead(i) w lines ls 4, \
     for [i=35:35]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=39:39]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=43:43]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=44:44]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=48:48]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=52:52]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=53:53]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=57:57]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=61:61]   filename using 1:(column(i))       title columnhead(i) w lines ls 5, \
     for [i=63:63]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=67:67]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=71:71]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=72:72]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=76:76]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=80:80]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=81:81]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=82:82]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=83:83]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=87:87]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=91:91]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=92:92]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=96:96]   filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=100:100] filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=101:101] filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=105:105] filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=109:109] filename using 1:(negsqrt(column(i))) title columnhead(i) w lines ls 1, \
     for [i=110:112] filename using 1:(column(i))       title columnhead(i) w lines ls 2
