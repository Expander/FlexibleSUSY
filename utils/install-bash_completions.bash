# This script installs bash completions for FlexibleSUSY's spectrum
# generators
#
# Usage:
#
#   . install-bash_completions.bash

_run_spectrum_generator()
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="
--slha-input-file=
--slha-output-file=
--spectrum-output-file=
--database-output-file=
--rgflow-output-file=
--build-info
--model-info
--help
--version
"

    # handle --xxxxxx=
    if [[ ${prev} == "--"* && ${cur} == "=" ]] ; then
        compopt -o filenames
        COMPREPLY=(*)
        return 0
    fi

    # handle --xxxxx=path
    if [[ ${prev} == '=' ]] ; then
        # unescape space
        cur=${cur//\\ / }
        # expand ~ to $HOME
        [[ ${cur} == "~/"* ]] && cur=${cur/\~/$HOME}
        # show completion if path exists (and escape spaces)
        compopt -o filenames
        local files=("${cur}"*)
        [[ -e ${files[0]} ]] && COMPREPLY=( "${files[@]// /\ }" )
        return 0
    fi

    # handle other options
    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    if [[ ${#COMPREPLY[@]} == 1 && ${COMPREPLY[0]} != "--"*"=" ]] ; then
        # if there's only one option, without =, then allow a space
        compopt +o nospace
    fi

    return 0
}

spectrum_generators=$(find models/ -name "run_*.x" -a ! -name "run_cmd_line*.x" -executable)

for sg in ${spectrum_generators}
do
    echo "installing bash completion for ${sg}"
    complete -o nospace -F _run_spectrum_generator "${sg}"
done
