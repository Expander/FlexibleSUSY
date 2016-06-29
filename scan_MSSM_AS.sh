start=0
stop=0.5
n_points=60

TB=5
Xt=0
MS=10000

./scan.sh --parameter=AS \
          --start=$start --stop=$stop --steps=$n_points \
          --AI="1.279440000e+04" --MZ=9.11876 \
          --TB="${TB}" --Xt="${Xt}" --MS="${MS}" \
    | tee scan_MSSM_AS.dat
