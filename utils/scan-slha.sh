#!/bin/sh

# This script allows one to scan over one parameter given an SLHA
# input file and a spectrum generator
#
# Examples:
#
# ./utils/scan-slha.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
#    --scan-range=MINPAR[1]=100-300:10 \
#    --output=MINPAR[1],MASS[25],Yu[3:3]
#
# cat model_files/CMSSM/LesHouches.in.CMSSM | \
# ./utils/scan-slha.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --scan-range=MINPAR[1]=100-300:10 \
#    --output=MINPAR[1],MASS[25],Yu[3:3]
#
# Author: Alexander Voigt

output=
scan_range=
slha_input=
slha_input_file=
spectrum_generator=

# prints SLHA block
print_slha_block_awk='
BEGIN {
   is_block = 0;
   if (block == "") {
      print "Error: block name not defined";
      print "   Please define the block name with -v block=<block-name>";
      exit 1
   }
}
{
   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$";
   not_pattern = "^block[[:blank:]]*.*$";

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   };

   if (is_block)
      print $0
}
'

# prints block entry
# expects block entry keys in the form x or x:y or x:y:z etc.
print_block_entry_awk='
{
  split(keys,k,":");

  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(length(k)+1)
}
'

# string to eval when script terminates, normally or not
actions_at_exit=""

#_____________________________________________________________________
at_exit() {
    local action
    for action in "$@"; do
	printf "%s\n" "$actions_at_exit" | grep -F -e "$action" > /dev/null ||
	actions_at_exit="$actions_at_exit
$action"
    done
}

#_____________________________________________________________________
do_actions_at_exit() {
    eval "$actions_at_exit"
}

#_____________________________________________________________________
help() {
cat <<EOF
Usage: ./`basename $0` [options]
Options:

  --output=             Comma separated list of output values
                        Syntax: <block>[<entry>]
                        Example: MINPAR[1],MASS[25],Yu[3:3]
  --scan-range=         Scan range
                        Syntax: <block>[<entry>]=<start>-<stop>:<steps-1>
                        Example: MINPAR[1]=100-1000:10
  --slha-input-file=    SLHA input file (optional).
                        If no SLHA input file is given, the SLHA input is
                        read from stdin .
  --spectrum-generator= Spectrum generator executable
  --help,-h             Print this help message

Examples:

   $ ./utils/scan-slha.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \\
        --scan-range=MINPAR[3]=1-30:21 \\
        --output=MINPAR[3],MASS[25],Yu[3:3] \\
     > scan-slha.dat

   $ cat model_files/CMSSM/LesHouches.in.CMSSM | \\
     ./utils/scan-slha.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --scan-range=MINPAR[3]=1-30:21 \\
        --output=MINPAR[3],MASS[25],Yu[3:3] \\
     > scan-slha.dat

   $  echo "set xlabel \"tan(beta)\"; 
            set ylabel \"mh / GeV\";
            plot 'scan-slha.dat' u 1:2 w linespoints" \\
        | gnuplot -p 
EOF
}

trap do_actions_at_exit 0
trap "exit 1" INT QUIT TERM

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --output=*)              output=$optarg ;;
            --scan-range=*)          scan_range=$optarg ;;
            --slha-input-file=*)     slha_input_file=$optarg ;;
            --spectrum-generator=*)  spectrum_generator=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

test -z "$slha_input_file" -o -e "$slha_input_file" || \
    { echo "Error: input file $slha_input_file not found." ; exit 1; }

slha_input=`cat ${slha_input_file}`

if test ! -x "$spectrum_generator"; then
    echo "Error: spectrum generator executable $spectrum_generator not found."
    exit 1
fi

if test -z "$output"; then
    echo "Error: no output values given.  Please provide them via --output= "
    exit 1
fi

if test -z "$scan_range"; then
    echo "Error: no scan range given.  Please provide it via --scan-range= "
    exit 1
fi

start=$(echo "$scan_range" | awk -F '[=:-]' '{ print $2 }')
stop=$(echo "$scan_range"  | awk -F '[=:-]' '{ print $3 }')
steps=$(echo "$scan_range" | awk -F : '{ print $NF }')
block=$(echo "$scan_range" | awk -F [ '{ print $1 }')
entry=$(echo "$scan_range" | awk -F '[][]' '{ print $2 }')

output_fields="$(echo $output | tr ',' ' ')"

# print comment line
if test "$steps" -gt 0; then
    comment="#"
    args=" "
    for f in $output_fields; do
        comment="$comment %16s"
        args="$args $f"
    done
    comment="$comment\n"
    printf "$comment" $args
fi

# start scan over points
for i in `seq 0 $steps`; do
    # calculate current value for the scanned variable
    value=$(cat <<EOF | bc
scale=20
$start + ($stop - $start)*${i} / $steps
EOF
    )

    # run the spectrum generator
    slha_output=$(
    { echo "$slha_input" ; \
      cat <<EOF
Block $block # added by `basename $0`
  $entry    $value
EOF
    } | $spectrum_generator --slha-input-file=- 2>/dev/null)

    printfstr=" "
    args=

    # get the output
    for f in $output_fields; do
        output_block=$(echo "$f" | awk -F [ '{ print $1 }')
        full_block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
        block_entries=$(echo "$f" | awk -F '[][]' '{ print $2 }')

        # get data value
        value=$(echo "$full_block" | awk -v keys="$block_entries" "$print_block_entry_awk")

        if test -z "$value"; then
            value="-"
        fi

        printfstr="$printfstr %16s"
        args="$args $value"
    done

    # print data line
    printfstr="$printfstr\n"
    printf "$printfstr" $args

done
