#!/bin/sh

[ $# -lt 5 ] && {
    echo "Usage: $0 <spectrum-generator> <input> <scan-fields> <output-fields> <mean-scale> [<factor>] [steps] [option]"
    echo "Option: --min --max --dif"
    exit 1
}

sg="$1"
shift
input="$1"
shift
scan_fields="$1"
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

print_what=dif

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --max) print_what=max ;;
            --min) print_what=min ;;
            --dif) print_what=dif ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

slha_input=$(cat ${input})

block=$(echo "$output_fields" | awk -F [ '{ print $1 }')
entry=$(echo "$output_fields" | awk -F '[][]' '{ print $2 }')
scan_block=$(echo "$scan_fields" | awk -F [ '{ print $1 }')
scan_entry=$(echo "$scan_fields" | awk -F '[][]' '{ print $2 }')

min_value_init=100000000
max_value_init=0
min_value="$min_value_init"
max_value="$max_value_init"

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
    scale=$(cat <<EOF | bc
scale = 15
start = ${mean_scale} / ${factor}
stop  = ${mean_scale} * ${factor}
start + (stop - start)*${i} / $steps
EOF
         )

    output=$({ echo "$slha_input";
      cat <<EOF
Block ${scan_block}
   ${scan_entry}  $scale   # renormalization scale
EOF
             } | "$sg" --slha-input-file=- 2> /dev/null)

    value=-

    [ -n "$output" ] && {
        value=$(echo "$output" | \
                       awk -f utils/print_slha_block.awk -v block="${block}" | \
                       awk -f utils/print_slha_block_entry.awk -v entries="${entry}" | tail -n 1)
    }

    [ -n "$value" -a "x$value" != "x-" ] && {
        min_value=$(min "$min_value" "$value" | sed -e 's/[eE]+*/*10^/')
        max_value=$(max "$max_value" "$value" | sed -e 's/[eE]+*/*10^/')
    }
done

[ "x$min_value" = "x$min_value_init" -o "x$max_value" = "x$max_value_init" ] && {
    min_value=0
    max_value=0
    echo "-"
    exit 1
}

delta=$(cat <<EOF | bc
scale = 15
define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

abs(${min_value} - ${max_value})
EOF
     )

case "$print_what" in
    max) echo "$max_value" ;;
    min) echo "$min_value" ;;
    dif) echo "$delta" ;;
    *)   echo "-" ;;
esac
