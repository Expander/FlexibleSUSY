#!/bin/sh

if [ $# -ne 4 ] ; then
    echo "Error: 4 Arguments required!"
    echo "  $0 <FH-dir> <MS> <tan(beta)> <Xt>"
    exit 1
fi

fh_dir="$1"
shift
MS="$1"
shift
TB="$1"
shift
Xt="$1"

fh="${fh_dir}/FeynHiggs"
fh_table="${fh_dir}/table"
fh_in=fh.in
fh_out="${fh_in}.fh-001"

At=$(echo "scale=16; $MS * $Xt + $MS / $TB" | bc -l)
Au=$(echo "scale=16; $MS/$TB" | bc -l)
Ab=$(echo "scale=16; $MS * $TB" | bc -l)

cat <<EOF > "${fh_in}"
Block MINPAR                 # Input parameters
    4   1.000000000e+00      # sign(mu)
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.384               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
BLOCK EXTPAR
         0     ${MS}   # Q
         1     ${MS}   # M1
         2     ${MS}   # M2
         3     ${MS}   # M3
        11     ${At}   # At
        12     ${Ab}   # Ab
        13     ${Ab}   # Atau
        23     ${MS}   # MUE
        26     ${MS}   # MA0
        25     ${TB}   # TB at Q
        31     ${MS}   # MSL(1)
        32     ${MS}   # MSL(2)
        33     ${MS}   # MSL(3)
        34     ${MS}   # MSE(1)
        35     ${MS}   # MSE(2)
        36     ${MS}   # MSE(3)
        41     ${MS}   # MSQ(1)
        42     ${MS}   # MSQ(2)
        43     ${MS}   # MSQ(3)
        44     ${MS}   # MSU(1)
        45     ${MS}   # MSU(2)
        46     ${MS}   # MSU(3)
        47     ${MS}   # MSD(1)
        48     ${MS}   # MSD(2)
        49     ${MS}   # MSD(3)
EOF

rm -f "${fh_out}"

${fh} "${fh_in}" 400203110 >/dev/null 2>&1

Mh=$(awk -f utils/print_slha_block.awk -v block=MASS "${fh_out}" | \
     awk -f utils/print_slha_block_entry.awk -v entries=25)

DMh=$(awk -f utils/print_slha_block.awk -v block=DMASS "${fh_out}" | \
      awk -f utils/print_slha_block_entry.awk -v entries=25)

[ "x$Mh" = "x" ] && Mh=-
[ "x$DMh" = "x" ] && DMh=-

echo "$Mh   $DMh"
