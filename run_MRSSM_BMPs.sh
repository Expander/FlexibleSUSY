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
    local block=
    local value=-
    local output_block=MASS
    local output_entry=25

    rm -f tmp.out tmp.in

    echo "$input" > tmp.in
    cat <<EOF >> tmp.in
Block SPhenoInput
     7   ${skip2L}            # Skip 2-loop Higgs corrections 
EOF

    ${sg} tmp.in tmp.out > /dev/null 2>&1

    if [ -f tmp.out ] ; then
        block=$(awk -v block="$output_block" "$print_slha_block_awk" tmp.out)
        value=$(echo "$block" | awk -v keys="$output_entry" "$print_block_entry_awk")
    fi

    rm -f tmp.in tmp.out

    echo "$value"
}

run_point_FlexibleSUSY() {
    local input="$1"
    local sg="$2"
    local value=-
    local output_block=MASS
    local output_entry=25

    slha_output=$(echo "$input" | ${sg} --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk")

    [ "x$value" = "x" ] && value="-"

    echo "$value"
}

max() {
    local M1="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local M2="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"

    echo $(cat <<EOF | bc -l
scale=10

define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

define max(i,j) {
    if (i > j) return i
    return j
}

max(abs($M1), abs($M2))
EOF
        )
}

calc_DMh_Q() {
    local Mh="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local MhQmax="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"
    local MhQmin="$(echo "$3" | sed -e 's/[eE]+*/*10^/')"

    echo $(cat <<EOF | bc -l
scale=10

define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

define min(i,j) {
    if (i < j) return i
    return j
}

define max(i,j) {
    if (i > j) return i
    return j
}

max(abs($Mh - $MhQmax), abs($Mh - $MhQmin))
EOF
        )
}

calc_DMh_4yt() {
    local M1="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local M2="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"
    local M3="$(echo "$3" | sed -e 's/[eE]+*/*10^/')"
    local M4="$(echo "$4" | sed -e 's/[eE]+*/*10^/')"

    echo $(cat <<EOF | bc -l
scale=10

define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

define max(i,j) {
    if (i > j) return i
    return j
}

define max6(a,b,c,d,e,f) {
    return max(a, max(b, max(c, max(d, max(e, f)))))
}

max6(abs($M1 - $M2), \
     abs($M1 - $M3), \
     abs($M1 - $M4), \
     abs($M2 - $M3), \
     abs($M2 - $M4), \
     abs($M3 - $M4)) / 2
EOF
        )
}

calc_quad_2() {
    local M1="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local M2="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"

    echo $(cat <<EOF | bc -l
scale=10
sqrt(($M1)^2 + ($M2)^2)
EOF
        )
}

absdiff() {
    local M1="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local M2="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"

    echo $(cat <<EOF | bc -l
define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

scale=10
abs($M1 - $M2)
EOF
        )
}

points="BM1p BM2p BM3p"

printf "#%6s %40s %40s %40s %40s %40s %40s\n" \
       "point" "SPheno-1L" "SPheno-2L" "SPheno-1L (FS-like)" "SPheno-2L (FS-like)" "FS-1L" "FStower-1L"

for p in ${points} ; do
    input=$(cat "LesHouches.in.MRSSM_${p}")

    input_SPheno_like="\
${input}
Block FlexibleSUSY
   17   1                  # Mt method (0 = FS, 1 = SPheno)
"

    input_yt_0L="\
${input}
Block FlexibleSUSY
   19   1    # mf tree-level matching
"

    MhSPheno1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2" "1")
    MhSPheno2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2" "0")
    MhSPhenoFS1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like" "1")
    MhSPhenoFS2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like" "0")
    # Yu at MZ FS-like
    MhFS1L=$(run_point_FlexibleSUSY "${input}" models/MRSSM2/run_MRSSM2.x)
    MhFStower1L=$(run_point_FlexibleSUSY "${input}" models/MRSSM2tower/run_MRSSM2tower.x)

    # Q uncertainty for SPheno
    DMhQmaxSPheno1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_uncertainty_max.sh" "1")
    DMhQminSPheno1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_uncertainty_min.sh" "1")
    DMhQmaxSPheno2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_uncertainty_max.sh" "0")
    DMhQminSPheno2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_uncertainty_min.sh" "0")

    # Q uncertainty for SPheno FS-like
    DMhQmaxSPhenoFS1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like_uncertainty_max.sh" "1")
    DMhQminSPhenoFS1L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like_uncertainty_min.sh" "1")
    DMhQmaxSPhenoFS2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like_uncertainty_max.sh" "0")
    DMhQminSPhenoFS2L=$(run_point_SPheno "${input}" "./SPhenoMRSSM2_FlexibleSUSY_like_uncertainty_min.sh" "0")

    # Q uncertainty for FS
    DMhQmaxFS1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2_uncertainty_max.sh")
    DMhQminFS1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2_uncertainty_min.sh")

    # Q uncertainty for FS-tower
    DMhQmaxFStower1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2tower_uncertainty_max.sh")
    DMhQminFStower1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2tower_uncertainty_min.sh")

    # Qmatch uncertainty for FS-tower
    DMhQmmaxFStower1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2tower_Qmatch_uncertainty_max.sh")
    DMhQmminFStower1L=$(run_point_FlexibleSUSY "${input}" "./MRSSM2tower_Qmatch_uncertainty_min.sh")

    # Yu at MZ SPheno-like
    MhFS1LSPhenoLike=$(run_point_FlexibleSUSY "${input_SPheno_like}" models/MRSSM2/run_MRSSM2.x)
    # Yu at MS FS-like
    MhFS1LYuatMS=$(run_point_FlexibleSUSY "${input}" models/MRSSM2YuatMS/run_MRSSM2YuatMS.x)
    # Yu at MS SPheno-like
    MhFS1LYuatMSSPhenoLike=$(run_point_FlexibleSUSY "${input_SPheno_like}" models/MRSSM2YuatMS/run_MRSSM2YuatMS.x)

    # tower, yt tree-level matching
    MhFStower1LYt0L=$(run_point_FlexibleSUSY "${input_yt_0L}" models/MRSSM2tower/run_MRSSM2tower.x)

    # DMh 4 x yt
    DMh4yt=$(calc_DMh_4yt "$MhFS1L" "$MhFS1LSPhenoLike" "$MhFS1LYuatMS" "$MhFS1LYuatMSSPhenoLike")

    # DMh SPheno Delta
    DMhSPheno1L=$(calc_quad_2 "$DMh4yt" \
                  "$(calc_DMh_Q "$MhSPheno1L" "$DMhQmaxSPheno1L" "$DMhQminSPheno1L")")

    DMhSPheno2L=$(calc_quad_2 "$DMh4yt" \
                  "$(calc_DMh_Q "$MhSPheno2L" "$DMhQmaxSPheno2L" "$DMhQminSPheno2L")")

    DMhSPhenoFS1L=$(calc_quad_2 "$DMh4yt" \
                  "$(calc_DMh_Q "$MhSPhenoFS1L" "$DMhQmaxSPhenoFS1L" "$DMhQminSPhenoFS1L")")

    DMhSPhenoFS2L=$(calc_quad_2 "$DMh4yt" \
                  "$(calc_DMh_Q "$MhSPhenoFS2L" "$DMhQmaxSPhenoFS2L" "$DMhQminSPhenoFS2L")")

    DMhFS1L=$(calc_quad_2 "$DMh4yt" \
             "$(calc_DMh_Q "$MhFS1L" "$DMhQmaxFS1L" "$DMhQminFS1L")")

    DMhFStower1L=$(calc_quad_2 \
                  $(max \
                      $(calc_DMh_Q "$MhFStower1L" "$DMhQmmaxFStower1L" "$DMhQmminFStower1L") \
                      $(absdiff    "$MhFStower1L" "$MhFStower1LYt0L")) \
                  $(calc_DMh_Q "$MhFStower1L" "$DMhQmaxFStower1L" "$DMhQminFStower1L"))

    printf " %6s %40s %40s %40s %40s %40s %40s\n" \
           "${p}" "${MhSPheno1L} +- ${DMhSPheno1L}" \
           "${MhSPheno2L} +- ${DMhSPheno2L}" \
           "${MhSPhenoFS1L} +- ${DMhSPhenoFS1L}" \
           "${MhSPhenoFS2L} +- ${DMhSPhenoFS2L}" \
           "${MhFS1L} +- ${DMhFS1L}" \
           "${MhFStower1L} +- ${DMhFStower1L}"
done
