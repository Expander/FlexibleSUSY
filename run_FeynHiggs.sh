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

At=$(echo "scale=16; $MS * $Xt + $MS / $TB" | bc -l)
Au=$(echo "scale=16; $MS/$TB" | bc -l)
Ab=$(echo "scale=16; $MS * $TB" | bc -l)

cat <<EOF > "${fh_in}"
invAlfaMZ 1.27944000E+02
AlfasMZ 1.184E-01
GF 1.16638E-05
MS 1.04E-1
MC 1.27
MT 173.34
MB 4.18
MZ 91.1876
TB $TB
MA0 $MS
MSusy $MS
Abs(MUE) $MS
Arg(MUE) 0
Abs(At) $At
Arg(At) 0
Abs(Ac) $Au
Arg(Ac) 0
Abs(Au) $Au
Arg(Au) 0
Abs(Ab) $Ab
Arg(Ab) 0
Abs(As) $Ab
Arg(As) 0
Abs(Ad) $Ab
Arg(Ad) 0
Abs(Atau) $Ab
Arg(Atau) 0
Abs(Amu) $Ab
Arg(Amu) 0
Abs(Ae) $Ab
Arg(Ae) 0
Abs(M_1) $MS
Arg(M_1) 0
Abs(M_2) $MS
Arg(M_2) 0
Abs(M_3) $MS
Arg(M_3) 0
scalefactor 1
EOF

${fh} "${fh_in}" 400203110 2>/dev/null | ${fh_table} Mh0 DeltaMh0 2>/dev/null
