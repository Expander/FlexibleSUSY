set terminal epslatex standalone size 12cm,8cm header '\renewcommand*{\familydefault}{\sfdefault}\usepackage[cm]{sfmath}\usepackage{amsmath}'

set xlabel "Energieskala / GeV"
set logscale x
set format x "$10^{%L}$"
set format y "$%g$"
set mytics 2
unset key
# set key box outside width 0.5

negsqrt(x) = x < 0 ? -sqrt(-x) : sqrt(x)

set style arrow 1 head filled size screen 0.02,10,45
set style arrow 1 head filled size screen 0.02,10,45

set style line 1 lt 1 lw 2 pt 7 ps 0.5 lc rgb "black"
set style line 2 lt 2 lw 2 pt 7 ps 0.5 lc rgb "red"
set style line 3 lt 3 lw 2 pt 7 ps 0.5 lc rgb "green"
set style line 4 lt 4 lw 2 pt 7 ps 0.5 lc rgb "blue"
set style line 5 lt 5 lw 2 pt 7 ps 0.5 lc rgb "magenta"

# parameter  column
# Mu         29
# BMu        62
# vu         34
# vd         35

##########

set output "MSSM_rgflow_0.tex"
filename='MSSM_rgflow_up_it0_constraint1'

set arrow 2 from 91,-1.2 to 91,-1 as 1 lt rgb "red" lw 2
set label 1 "$\\color{red}{M_Z}$"        at 91,-1.3 center
set label 6 "$g_1$" at 20,0.457664 center
set label 7 "$g_2$" at 20,0.641152 center
set label 8 "$g_3$" at 20,1.11384  center

plot [1:1.0e18] [-1:1.25] \
     '< head -n 2 '.filename using 1:30 title "$g_i$" with points ls 1, \
     '< head -n 2 '.filename using 1:31 title ""      with points ls 1, \
     '< head -n 2 '.filename using 1:32 title ""      with points ls 1

unset arrow 2
unset label 1

set output
system 'pdflatex MSSM_rgflow_0.tex'

##########

set output "MSSM_rgflow_1.tex"
filename='MSSM_rgflow_up_it0_constraint1'

set arrow 1 from 1.e8,1     to 1.e10,1  as 1
set arrow 2 from 2e+16,-1.2 to 2e+16,-1 as 1 lt rgb "red" lw 2
set label 1 "$\\color{red}{M_X}$"       at 2e+16,-1.3 center
set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73

plot [1:1.0e18] [-1:1.25] \
     filename using 1:30 title "$g_i$" with lines ls 1, \
     filename using 1:31 title ""      with lines ls 1, \
     filename using 1:32 title ""      with lines ls 1

unset arrow 1
unset arrow 2
unset label 1
unset label 2
unset label 6
unset label 7
unset label 8

set output
system 'pdflatex MSSM_rgflow_1.tex'

##########

set arrow 2 from 2e+16,-1.2 to 2e+16,-1 as 1 lt rgb "red" lw 2
set label 1 "$\\color{red}{M_X}$"       at 2e+16,-1.3 center
set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_2.tex"
filename='MSSM_rgflow_down_it0_constraint1'

plot [1:1.0e18] [-1:1.25] \
     'MSSM_rgflow_up_it0_constraint1' using 1:30 title "$g_i$"       w lines ls 1, \
     'MSSM_rgflow_up_it0_constraint1' using 1:31 title ""            w lines ls 1, \
     'MSSM_rgflow_up_it0_constraint1' using 1:32 title ""            w lines ls 1, \
     for [i=35:35]   '< tail -n 1 '.filename using 1:(column(i)/1000)          title "$T_i$" w points ls 2, \
     for [i=63:63]   '< tail -n 1 '.filename using 1:(negsqrt(column(i))/1000) title "$m_i$" w points ls 3, \
     for [i=110:110] '< tail -n 1 '.filename using 1:(column(i)/1000)          title "$M_i$" w points ls 4

unset arrow 2
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_2.tex'

##########

set arrow 1 from 1.e10,1 to 1.e8,1 as 1
set arrow 2 from 1007.78,-1.2 to 1007.78,-1 as 1 lt rgb "red" lw 2
set label 1 "$\\color{red}{M_S}$" at 1007.78,-1.3 center

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_3.tex"
filename='MSSM_rgflow_down_it0_constraint1'

plot [1:1.0e18] [-1:1.25] \
                     filename using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename using 1:31                        title ""            w lines ls 1, \
                     filename using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename using 1:(column(i)/1000)          title ""            w lines ls 4

unset arrow 1
unset arrow 2
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_3.tex'

##########

set arrow 1 from 1.e10,1 to 1.e8,1 as 1
set arrow 2 from 91.1876,-1.2 to 91.1876,-1 as 1 lt rgb "red" lw 2
set label 1 "$\\color{red}{M_\\text{Z}}$" at 91.1876,-1.3 center

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_4.tex"
filename2='MSSM_rgflow_down_it0_constraint2'

plot [1:1.0e18] [-1:1.25] \
                     filename  using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename  using 1:31                        title ""            w lines ls 1, \
                     filename  using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename  using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename  using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename  using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename  using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename  using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename  using 1:(column(i)/1000)          title ""            w lines ls 4, \
                     filename2 using 1:30                        title ""            w lines ls 1, \
                     filename2 using 1:31                        title ""            w lines ls 1, \
                     filename2 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=39:39]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=67:67]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:112] filename2 using 1:(column(i)/1000)          title ""            w lines ls 4, \
     for [i=29:29]   filename2 using 1:(column(i)/1000)          title "$\\mu$"      w lines ls 5, \
     for [i=62:62]   filename2 using 1:(negsqrt(column(i))/1000) title "$B\\mu$"     w lines ls 5

unset arrow 1
unset arrow 2
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_4.tex'

########## iteration 1

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_5.tex"
filename1='MSSM_rgflow_down_it1_constraint1'
filename2='MSSM_rgflow_down_it1_constraint2'

plot [1:1.0e18] [-1:1.25] \
                     filename1 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename1 using 1:31                        title ""            w lines ls 1, \
                     filename1 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename1 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename1 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename1 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename1 using 1:(column(i)/1000)          title ""            w lines ls 4, \
                     filename2 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename2 using 1:31                        title ""            w lines ls 1, \
                     filename2 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename2 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename2 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename2 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename2 using 1:(column(i)/1000)          title ""            w lines ls 4, \
     for [i=29:29]   filename2 using 1:(column(i)/1000)          title "$\\mu$"      w lines ls 5, \
     for [i=62:62]   filename2 using 1:(negsqrt(column(i))/1000) title "$B\\mu$"     w lines ls 5

unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_5.tex'

########## iteration 2

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_6.tex"
filename1='MSSM_rgflow_down_it2_constraint1'
filename2='MSSM_rgflow_down_it2_constraint2'

plot [1:1.0e18] [-1:1.25] \
                     filename1 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename1 using 1:31                        title ""            w lines ls 1, \
                     filename1 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename1 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename1 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename1 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename1 using 1:(column(i)/1000)          title ""            w lines ls 4, \
                     filename2 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename2 using 1:31                        title ""            w lines ls 1, \
                     filename2 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename2 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename2 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename2 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename2 using 1:(column(i)/1000)          title ""            w lines ls 4, \
     for [i=29:29]   filename2 using 1:(column(i)/1000)          title "$\\mu$"      w lines ls 5, \
     for [i=62:62]   filename2 using 1:(negsqrt(column(i))/1000) title "$B\\mu$"     w lines ls 5

unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_6.tex'

########## iteration 7

set label 2 "\\scalebox{0.7}{$g_0$}"     at 3.e16,0.73
set label 3 "\\scalebox{0.7}{$m_0$}"     at 3.e16,0.125
set label 4 "\\scalebox{0.7}{$M_{1/2}$}" at 3.e16,0.5
set label 5 "\\scalebox{0.7}{$A_0$}"     at 3.e16,0

set output "MSSM_rgflow_7.tex"
filename1='MSSM_rgflow_down_it25_constraint1'
filename2='MSSM_rgflow_down_it25_constraint2'

plot [1:1.0e18] [-1:1.25] \
                     filename1 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename1 using 1:31                        title ""            w lines ls 1, \
                     filename1 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename1 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename1 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename1 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename1 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename1 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename1 using 1:(column(i)/1000)          title ""            w lines ls 4, \
                     filename2 using 1:30                        title "$g_i$"       w lines ls 1, \
                     filename2 using 1:31                        title ""            w lines ls 1, \
                     filename2 using 1:32                        title ""            w lines ls 1, \
     for [i=35:35]   filename2 using 1:(column(i)/1000)          title "$T_i$"       w lines ls 2, \
     for [i=39:39]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=43:43]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=44:44]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=48:48]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=52:52]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=53:53]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=57:57]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=61:61]   filename2 using 1:(column(i)/1000)          title ""            w lines ls 2, \
     for [i=63:63]   filename2 using 1:(negsqrt(column(i))/1000) title "$m_i$"       w lines ls 3, \
     for [i=67:67]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=71:71]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=72:72]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=76:76]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=80:80]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=81:81]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=82:82]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=83:83]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=87:87]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=91:91]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=92:92]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=96:96]   filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=100:100] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=101:101] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=105:105] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=109:109] filename2 using 1:(negsqrt(column(i))/1000) title ""            w lines ls 3, \
     for [i=110:110] filename2 using 1:(column(i)/1000)          title "$M_i$"       w lines ls 4, \
     for [i=111:112] filename2 using 1:(column(i)/1000)          title ""            w lines ls 4, \
     for [i=29:29]   filename2 using 1:(column(i)/1000)          title "$\\mu$"      w lines ls 5, \
     for [i=62:62]   filename2 using 1:(negsqrt(column(i))/1000) title "$B\\mu$"     w lines ls 5

unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set output
system 'pdflatex MSSM_rgflow_7.tex'
