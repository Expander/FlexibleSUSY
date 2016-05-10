#!/bin/sh

[ $# -lt 4 ] && {
    echo "Usage: $0 <SPheno> <input-file> <output-fields> <mean-scale> [<factor>] [steps]"
    exit 1
}

sg="$1"
shift
input="$1"
shift
output_fields="$1"
shift
mean_scale=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
shift

factor=2
[ "$#" -ge 1 ] && {
    factor="$1"
    shift
}

steps=10
[ "$#" -ge 1 ] && {
    steps="$1"
    shift
}

input_tmp="${input}.$$"
output="SPheno.out.$$"
block=$(echo "$output_fields" | awk -F [ '{ print $1 }')
entry=$(echo "$output_fields" | awk -F '[][]' '{ print $2 }')

min_value=100000000
max_value=0

# returns minimum of two numbers
min() {
    local a=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    local b=$(echo "$2" | sed -e 's/[eE]+*/*10^/')
    cat <<EOF | bc -l
define min(i,j) {
    if (i < j) return i
    return j
}

min($a,$b)
EOF
}

# returns maximum of two numbers
max() {
    local a=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    local b=$(echo "$2" | sed -e 's/[eE]+*/*10^/')
    cat <<EOF | bc -l
define max(i,j) {
    if (i > j) return i
    return j
}

max($a,$b)
EOF
}

for i in $(seq 0 ${steps}); do
    rm -f "$output" "$input_tmp"

    scale=$(cat <<EOF | bc
scale = 15
start = ${mean_scale} / ${factor}
stop  = ${mean_scale} * ${factor}
start + (stop - start)*${i} / $steps
EOF
         )

    { cat "$input";
      cat <<EOF
Block SPhenoInput
   33  $scale   # renormalization scale
EOF
    } > "$input_tmp"

    "$sg" "$input_tmp" "$output" > /dev/null

    value=-

    [ -f "$output" ] && {
        value=$(awk -f utils/print_slha_block.awk -v block="${block}" "$output" | \
                awk -f utils/print_slha_block_entry.awk -v entries="${entry}" | tail -n 1)
    }

    [ "x$value" != "x-" ] && {
        min_value=$(min "$min_value" "$value" | sed -e 's/[eE]+*/*10^/')
        max_value=$(max "$max_value" "$value" | sed -e 's/[eE]+*/*10^/')
    }
done

rm -f "$output" "$input_tmp"

delta=$(cat <<EOF | bc
scale = 15
define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

abs(${min_value} - ${max_value})
EOF
        )

echo "$delta"
