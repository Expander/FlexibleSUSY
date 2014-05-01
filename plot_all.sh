
# MSSM

# root -l -b -q higgs-study/plot2d.C+(\"higgs-study/data/scan_MSSM.dat\",\"CMSSM \$m_{h}$ / GeV\")

# NMSSM

# root -l -b -q higgs-study/plot2d.C+(\"higgs-study/data/scan_NMSSM_lambda-0.1.dat\",\"CNMSSM \$m_{h}$ / GeV\")

# root -l -b -q higgs-study/plot2d.C+(\"higgs-study/data/scan_NMSSM_lambda-0.2.dat\",\"CNMSSM \$m_{h}$ / GeV\")

files="`ls higgs-study/data/scan_*.dat`"

for f in ${files}
do
    root -l -b -q higgs-study/plot2d.C+(\"${f}\",\"\$m_{h}$ / GeV\")
done
