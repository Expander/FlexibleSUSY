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
  len = split(keys,k,":");

  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(len + 1)
}
'

run_point_SPheno() {
    local input="$1"
    local sg="$2"
    local skip2L="$3"
    local value=

    rm -f tmp.out tmp.in

    cp "$input" tmp.in
    cat <<EOF >> tmp.in
Block SPhenoInput
     7   ${skip2L}            # Skip 2-loop Higgs corrections 
EOF

    ${sg} tmp.in tmp.out > /dev/null 2>&1
    [ -f tmp.out ] && value=$(grep hh_1 tmp.out | awk '{print $2}')

    rm -f tmp.in tmp.out

    echo "$value"
}

run_point_FlexibleSUSY() {
    local input="$1"
    local sg="$2"
    local value=-
    local output_block=MASS
    local output_entry=25

    slha_output=$(cat "$input" | ${sg} --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk")

    [ "x$value" = "x" ] && value="-"

    echo "$value"
}

points="BM1p BM2p BM3p"

printf "# %20s %20s %20s %20s %20s %20s %20s\n" "point" "SPheno-1L" "SPheno-2L" "SPheno-1L (FS-like)" "SPheno-2L (FS-like)" "FS-1L" "FStower-1L"

for p in ${points} ; do
    input_file="LesHouches.in.MRSSM_${p}"

    MhSPheno1L=$(run_point_SPheno "${input_file}" "./SPhenoMRSSM2" "1")
    MhSPheno2L=$(run_point_SPheno "${input_file}" "./SPhenoMRSSM2" "0")
    MhSPhenoFS1L=$(run_point_SPheno "${input_file}" "./SPhenoMRSSM2_FlexibleSUSY_like" "1")
    MhSPhenoFS2L=$(run_point_SPheno "${input_file}" "./SPhenoMRSSM2_FlexibleSUSY_like" "0")
    MhFS1L=$(run_point_FlexibleSUSY "${input_file}" models/MRSSM2/run_MRSSM2.x)
    MhFStower1L=$(run_point_FlexibleSUSY "${input_file}" models/MRSSM2tower/run_MRSSM2tower.x)

    printf "  %20s %20s %20s %20s %20s %20s %20s\n" "${p}" "${MhSPheno1L}" "${MhSPheno2L}"  "${MhSPhenoFS1L}" "${MhSPhenoFS2L}" "${MhFS1L}" "${MhFStower1L}"
done
