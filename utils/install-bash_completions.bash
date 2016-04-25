# This script installs bash completions for FlexibleSUSY's configure,
# createmodel, createaddon scripts as well as all spectrum generators.
#
# Usage:
#
#   . install-bash_completions.bash

__find_available_models()
{
    find $(dirname $0)/model_files/ -mindepth 2 -type f -name FlexibleSUSY.m.in -exec bash -c 'basename $(dirname $1)' modelname {} \; | sort -u
}

__find_created_models()
{
    find $(dirname $0)/models/ -mindepth 2 -type f -name module.mk -exec bash -c 'basename $(dirname $1)' modelname {} \; | sort -u
}

__find_available_addons()
{
    find $(dirname $0)/addons/ -mindepth 1 -type d -name "*" -exec basename {} \; | sort -u
}

__build_filename_completion_list()
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="$@"

    [ "x$opts" = x ] && return 0

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

__select_from_list()
{
    local pprev="$1"
    shift
    local prev="$1"
    shift
    local cur="$1"
    shift
    local option="$1"
    shift
    local list="$@"

    # handle --XXX=
    if [[ ${prev}${cur} == "${option}" ]] ; then
        echo "$list"
        return 0
    fi

    # handle --XXX=YYY
    if [[ ${pprev}${prev} == "${option}" ]] ; then
        local x
        local selected_list=$(for x in ${list} ; do [[ $x == "${cur}"* ]] && echo $x; done)
        echo "$selected_list"
        return 0
    fi
}

_run_spectrum_generator()
{
    local opts="
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

    __build_filename_completion_list "${opts}"
}

_configure()
{
    local cur prev pprev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    pprev="${COMP_WORDS[COMP_CWORD-2]}"
    opts="
--enable-colors
--enable-compile
--enable-compiler-warnings
--enable-debug
--enable-mass-error-check
--enable-fflite
--enable-looptools
--enable-meta
--enable-silent
--enable-sqlite
--enable-static-libs
--enable-threads
--enable-tsil
--enable-verbose

--disable-colors
--disable-compile
--disable-compiler-warnings
--disable-debug
--disable-mass-error-check
--disable-fflite
--disable-looptools
--disable-meta
--disable-silent
--disable-sqlite
--disable-static-libs
--disable-threads
--disable-tsil
--disable-verbose

--with-addons=
--with-algorithms=
--with-blas-libdir=
--with-blas-libs=
--with-boost-libdir=
--with-boost-incdir=
--with-cxx=
--with-cxxflags=
--with-eigen-incdir=
--with-fc=
--with-fflags=
--with-flibs=
--with-gsl-config=
--with-install-dir=
--with-lapack-libdir=
--with-lapack-libs=
--with-ldflags=
--with-ldlibs=
--with-lib-ext=
--with-looptools-libdir=
--with-looptools-incdir=
--with-make-lib-cmd=
--with-math-cmd=
--with-models=
--with-optional-modules=
--with-sqlite-libdir=
--with-sqlite-incdir=
--with-tsil-libdir=
--with-tsil-incdir=

--help
--version
"

    # handle --with-addons=
    if [[ ${prev}${cur} == "--with-addons=" || ${pprev}${prev} == "--with-addons=" ]] ; then
        local available_addons=$(__find_available_addons)
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-addons=" "$available_addons") )
        return 0
    fi

    # handle --with-algorithms=
    if [[ ${prev}${cur} == "--with-algorithms=" || ${pprev}${prev} == "--with-algorithms=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-algorithms=" "two_scale lattice") )
        return 0
    fi

    # handle --with-cxx=
    if [[ ${prev}${cur} == "--with-cxx=" || ${pprev}${prev} == "--with-cxx=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-cxx=" $(compgen -c)) )
        return 0
    fi

    # handle --with-cxxflags=
    if [[ ${prev}${cur} == "--with-cxxflags=" || ${pprev}${prev} == "--with-cxxflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-fc=
    if [[ ${prev}${cur} == "--with-fc=" || ${pprev}${prev} == "--with-fc=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-fc=" $(compgen -c)) )
        return 0
    fi

    # handle --with-fflags=
    if [[ ${prev}${cur} == "--with-fflags=" || ${pprev}${prev} == "--with-fflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-gsl-config=
    if [[ ${prev}${cur} == "--with-gsl-config=" || ${pprev}${prev} == "--with-gsl-config=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-gsl-config=" $(compgen -c)) )
        return 0
    fi

    # handle --with-ldflags=
    if [[ ${prev}${cur} == "--with-ldflags=" || ${pprev}${prev} == "--with-ldflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-libext=
    if [[ ${prev}${cur} == "--with-libext=" || ${pprev}${prev} == "--with-libext=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-make-lib-cmd=
    if [[ ${prev}${cur} == "--with-make-lib-cmd=" || ${pprev}${prev} == "--with-make-lib-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-make-lib-cmd=" $(compgen -c)) )
        return 0
    fi

    # handle --with-math-cmd=
    if [[ ${prev}${cur} == "--with-math-cmd=" || ${pprev}${prev} == "--with-math-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-math-cmd=" $(compgen -c)) )
        return 0
    fi

    # handle --with-models=
    if [[ ${prev}${cur} == "--with-models=" || ${pprev}${prev} == "--with-models=" ]] ; then
        local available_models=$(__find_created_models)
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-models=" "$available_models") )
        return 0
    fi

    # handle --with-optional-modules=
    if [[ ${prev}${cur} == "--with-optional-modules=" || ${pprev}${prev} == "--with-optional-modules=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-optional-modules=" "examples test") )
        return 0
    fi

    __build_filename_completion_list "${opts}"
}

_createaddon()
{
    local cur opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    opts="
--name=
--force
--help
"

    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    if [[ ${#COMPREPLY[@]} == 1 && ${COMPREPLY[0]} != "--"*"=" ]] ; then
        # if there's only one option, without =, then allow a space
        compopt +o nospace
    fi

    return 0
}

_createmodel()
{
    local cur prev pprev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    pprev="${COMP_WORDS[COMP_CWORD-2]}"
    opts="
--name=
--model-file=
--sarah-model=
--with-math-cmd=
--enable-debug
--force
--help
"

    # handle --model-file=
    if [[ ${prev}${cur} == "--model-file=" || ${pprev}${prev} == "--model-file=" ]] ; then
        local available_models=$(__find_available_models)
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--model-file=" "$available_models") )
        return 0
    fi

    # handle --with-math-cmd=
    if [[ ${prev}${cur} == "--with-math-cmd=" || ${pprev}${prev} == "--with-math-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-math-cmd=" $(compgen -c)) )
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

echo "installing bash completion for configure"
complete -o nospace -F _configure configure

echo "installing bash completion for createmodel"
complete -o nospace -F _createmodel createmodel

echo "installing bash completion for createaddon"
complete -o nospace -F _createaddon createaddon

spectrum_generators=$(find models/ -name "run_*.x" -a ! -name "run_cmd_line*.x" -executable)

for sg in ${spectrum_generators}
do
    echo "installing bash completion for ${sg}"
    complete -o nospace -F _run_spectrum_generator "${sg}"
done
