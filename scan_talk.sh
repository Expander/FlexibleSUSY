./scan.sh --parameter=MS --start=91 --stop=1000   --steps=100 --step-size=linear --TB=5  --Xt=0     > ~/talks/KUTS-2015-Heidelberg/plots/scale_low_TB5.dat
./scan.sh --parameter=MS --start=91 --stop=1000   --steps=100 --step-size=linear --TB=20 --Xt=0     > ~/talks/KUTS-2015-Heidelberg/plots/scale_low_TB20.dat
./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log    --TB=5  --Xt=0     > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB5.dat
./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log    --TB=20 --Xt=0     > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB20.dat
./scan.sh --parameter=Xt --start=-3 --stop=3      --steps=60  --step-size=linear --TB=5  --MS=1000  > ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS1000.dat
./scan.sh --parameter=Xt --start=-3 --stop=3      --steps=60  --step-size=linear --TB=5  --MS=10000 > ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS10000.dat
./scan.sh --parameter=Xt --start=-3 --stop=3      --steps=60  --step-size=linear --TB=20 --MS=2000  > ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB20_MS2000.dat

./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=5  --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB5_Xtsqrt6.dat
./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=20 --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB20_Xtsqrt6.dat

# with tree-level matching + 1L lambda matching
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=100   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS100_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=200   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS200_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=300   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS300_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=400   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS400_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=500   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS500_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=1000  | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS1000_Matching0L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=10000 | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS10000_Matching0L.dat

./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=5  --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB5_Xtsqrt6_Matching0L.dat
./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=20 --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB20_Xtsqrt6_Matching0L.dat

# with tree-level matching + 1L lambda matching + 1L pole masses

./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=100   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS100_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=200   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS200_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=300   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS300_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=400   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS400_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=500   | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS500_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=1000  | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS1000_Matching0L_SE1L.dat
./scan.sh --parameter=Xt --start=-3 --stop=3   --steps=100 --step-size=linear --TB=5 --MS=10000 | tee ~/talks/KUTS-2015-Heidelberg/plots/Xt_TB5_MS10000_Matching0L_SE1L.dat

./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=5  --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB5_Xtsqrt6_Matching0L_SE1L.dat
./scan.sh --parameter=MS --start=91 --stop=100000 --steps=100 --step-size=log --TB=20 --Xt=2.44949 > ~/talks/KUTS-2015-Heidelberg/plots/scale_high_TB20_Xtsqrt6_Matching0L_SE1L.dat
