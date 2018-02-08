#!/bin/sh

# This script provides the get_model_specific_corrections() function,
# which takes a FlexibleSUSY model file [$1] (and optionally the
# Mathematica kernel command [$2]) and returns the list of
# model-specific corrections necessary for the build.

#_____________________________________________________________________
__contains() {
    # Check if string $1 contains $2
    if test $# -lt 2 ; then
        echo "contains: Too few arguments"
        return 1
    fi
    local string="$1"
    local substring="$2"
    for f in ${string}; do
        test "x$f" = "x$substring" && return 0
    done
    return 1 # not found
}

#_____________________________________________________________________
__make_unique() {
    # removes duplicate words from the string $1 (preserves order)
    local str="$1"
    local res=""

    for s in ${str}; do
        __contains "${res}" "${s}" && continue
        res="${res} ${s}"
    done

    echo "${res}" | sed -e 's/^ *//' -e 's/ *$//'
}

#_____________________________________________________________________
__mma_symbol_is() {
    local model_file="$1"
    local math_cmd="$2"
    local symbol="$3"
    local value="$4"

    [ -e "${model_file}" ] || {
        >&2 echo "Error: file \"${model_file}\" does not exist!"
        return 1
    }

    command -v "${math_cmd}" >/dev/null 2>&1 || {
        >&2 echo "Error: Mathematica kernel \"${math_cmd}\" is not executable!"
        return 1
    }

    [ -z "${symbol}" ] && {
        >&2 echo "Error: no symbol given!"
        return 1
    }

    [ -z "${value}" ] && {
        >&2 echo "Error: no value given!"
        return 1
    }

    echo "Get[\"${model_file}\"]; Quit[Boole[!(${symbol} === ${value})]]" | "${math_cmd}" > /dev/null 2>&1

    return $?
}

#_____________________________________________________________________
get_model_specific_corrections() {
    local model_file="$1"
    local math_cmd="${2:-math}"

    # always build SM, because it is needed for the effective couplings module
    local ho="SM"

    [ -e "${model_file}" ] || {
        >&2 echo "Error: file \"${model_file}\" does not exist!"
        echo "${ho}"
        return 1
    }

    command -v "${math_cmd}" >/dev/null 2>&1 || {
        >&2 echo "Error: Mathematica kernel \"${math_cmd}\" is not executable!"
        >&2 echo "   Please specify the Mathematica kernel as 2nd argument!"
        echo "${ho}"
        return 1
    }

    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs2LoopSM" "True"    && ho="${ho} SM"
    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs3LoopSM" "True"    && ho="${ho} SM"
    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs4LoopSM" "True"    && ho="${ho} SM"
    __mma_symbol_is "$model_file" "$math_cmd" "UseSMAlphaS3Loop" "True"   && ho="${ho} SM"
    __mma_symbol_is "$model_file" "$math_cmd" "FlexibleEFTHiggs" "True"   && ho="${ho} SM"

    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs2LoopMSSM" "True"  && ho="${ho} MSSM_higgs"
    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs3LoopMSSM" "True"  && ho="${ho} MSSM_higgs"
    __mma_symbol_is "$model_file" "$math_cmd" "UseMSSMAlphaS2Loop" "True" && ho="${ho} MSSM_thresholds"
    __mma_symbol_is "$model_file" "$math_cmd" "UseMSSMYukawa2Loop" "True" && ho="${ho} MSSM_thresholds"

    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs2LoopNMSSM" "True" && ho="${ho} MSSM_higgs NMSSM_higgs"

    __mma_symbol_is "$model_file" "$math_cmd" "UseHiggs3LoopSplit" "True" && ho="${ho} SplitMSSM"

    echo "$(__make_unique "${ho}")"
}
