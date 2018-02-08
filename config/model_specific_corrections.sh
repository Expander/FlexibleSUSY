#!/bin/sh

#_____________________________________________________________________
__make_unique() {
    # removes duplicate words from the string $1 (preserves order)
    local str="$1"
    local res=""

    for s in ${str}; do
        contains "${res}" "${s}" && continue
        res="${res} ${s}"
    done

    echo "${res}" | sed -e 's/^ *//' -e 's/ *$//'
}

#_____________________________________________________________________
get_model_specific_corrections() {
    local model_file="$1"

    # always build SM, because it is needed for the effective couplings module
    local ho="SM"

    if test -e "$model_file"; then
        grep "UseHiggs2LoopSM *= *True"    "$model_file" > /dev/null 2>&1 && ho="${ho} SM"
        grep "UseHiggs3LoopSM *= *True"    "$model_file" > /dev/null 2>&1 && ho="${ho} SM"
        grep "UseHiggs4LoopSM *= *True"    "$model_file" > /dev/null 2>&1 && ho="${ho} SM"
        grep "UseSMAlphaS3Loop *= *True"   "$model_file" > /dev/null 2>&1 && ho="${ho} SM"
        grep "FlexibleEFTHiggs *= *True"   "$model_file" > /dev/null 2>&1 && ho="${ho} SM"

        grep "UseHiggs2LoopMSSM *= *True"  "$model_file" > /dev/null 2>&1 && ho="${ho} MSSM_higgs"
        grep "UseHiggs3LoopMSSM *= *True"  "$model_file" > /dev/null 2>&1 && ho="${ho} MSSM_higgs"
        grep "UseMSSMAlphaS2Loop *= *True" "$model_file" > /dev/null 2>&1 && ho="${ho} MSSM_thresholds"
        grep "UseMSSMYukawa2Loop *= *True" "$model_file" > /dev/null 2>&1 && ho="${ho} MSSM_thresholds"

        grep "UseHiggs2LoopNMSSM *= *True" "$model_file" > /dev/null 2>&1 && ho="${ho} MSSM_higgs NMSSM_higgs"

        grep "UseHiggs3LoopSplit *= *True" "$model_file" > /dev/null 2>&1 && ho="${ho} SplitMSSM"
    fi

    echo "$(__make_unique "${ho}")"
}
