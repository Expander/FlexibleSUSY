#!/bin/sh

[ $# -lt 1 ] && {
    echo "Usage: $0 --slha-input-file=<input-file>"
    exit 1
}

input=-
output_block=MASS
output_entry=25
scale_block=EXTPAR
scale_entry=103

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac
        case $1 in
            --slha-input-file=*)     input=$optarg ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

slha_input=$(cat ${input})

mean_scale=$(echo "$slha_input" | \
                    awk -f utils/print_slha_block.awk -v block="MINPAR" | \
                    awk -f utils/print_slha_block_entry.awk -v entries="1" | tail -n 1)

[ -z "$mean_scale" ] && {
    echo "Error: mean scale not found in MINPAR[1] in the following input:" > /dev/stderr
    echo "$slha_input" > /dev/stderr
    exit 1
}

delta=$(echo "$slha_input" | ./scan_scale_FlexibleSUSY.sh models/MRSSM2tower/run_MRSSM2tower.x - "${scale_block}[${scale_entry}]" "${output_block}[${output_entry}]" "$mean_scale" 2 10 --max)

echo "$slha_input"
cat <<EOF
Block ${output_block}
   ${output_entry}  ${delta}
EOF
