
files="`ls higgs-study/data/scan_*.dat`"

for f in ${files}
do
    root -l -b -q higgs-study/plot2d.C+(\"${f}\",\"\$m_{h}$ / GeV\")
done
