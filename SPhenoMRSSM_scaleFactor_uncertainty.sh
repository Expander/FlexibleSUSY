#!/bin/sh

[ $# -lt 2 ] && {
    echo "Usage: $0 <input-file> <output-file>"
    exit 1
}

input="$1"
shift
output="$1"

output_block=MASS
output_entry=25
scale_block=MINPAR
scale_entry=1

mean_scale=$(awk -f utils/print_slha_block.awk -v block="$scale_block" "$input" | \
             awk -f utils/print_slha_block_entry.awk -v entries="$scale_entry" | tail -n 1)

delta=$(./scan_scale_SPheno.sh ./SPhenoMRSSM_scaleFactor "$input" "${output_block}[${output_entry}]" "$mean_scale")

block=$(cat <<EOF
Block ${output_block}
   ${output_entry}  ${delta}
EOF
     )

if [ -f "${output}" ] ; then
    echo "$block" >> "$output"
else
    echo "$block" > "$output"
fi
