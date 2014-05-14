set terminal x11

set title ""
set xlabel '$t_c$'
# set xlabel '$\tan\beta$'
set ylabel '$m_h^{\text{pole}}$'

set key box bottom right

if (!exists("filename")) filename='higgs-study/data/tc/scan_MSSM_tc.dat'

# plot "<awk '{ if ($2 == 0    && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0$'   , \
#      "<awk '{ if ($2 == 0.01 && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0.01$', \
#      "<awk '{ if ($2 == 0.02 && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0.02$', \
#      "<awk '{ if ($2 == 0.03 && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0.03$', \
#      "<awk '{ if ($2 == 0.04 && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0.04$', \
#      "<awk '{ if ($2 == 0.05 && $4 == 0) print $0 }' ".filename using 1:3 title '$t_c = 0.05$'

plot "<awk '{ if ($1 == 5  && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 5$' , \
     "<awk '{ if ($1 == 10 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 10$', \
     "<awk '{ if ($1 == 15 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 15$', \
     "<awk '{ if ($1 == 20 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 20$', \
     "<awk '{ if ($1 == 25 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 25$', \
     "<awk '{ if ($1 == 30 && $4 == 0) print $0 }' ".filename using 2:3 title '$\tan\beta = 30$'
