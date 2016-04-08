# This script installs bash completions for FlexibleSUSY's configure,
# createmodel, createaddon scripts as well as all spectrum generators.
#
# Usage:
#
#   . install-bash_completions.bash

_build_completion_list()
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

    _build_completion_list "${opts}"
}

_configure()
{
    local opts="
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
--with-optional-modules=
--with-sqlite-libdir=
--with-sqlite-incdir=
--with-tsil-libdir=
--with-tsil-incdir=
--with-math-cmd=
--with-models=

--help
--version
"

    _build_completion_list "${opts}"
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

    local available_models=$(find $(dirname $0)/model_files/ -mindepth 2 -type f -name FlexibleSUSY.m.in -exec bash -c 'basename $(dirname $1)' modelname {} \; | sort -u)

    # handle --model-file=
    if [[ ${prev} == "--model-file" && ${cur} == "=" ]] ; then
        COMPREPLY=( $available_models )
        return 0
    fi

    # handle --model-file=XXX
    if [[ ${pprev} == "--model-file" && ${prev} == "=" ]] ; then
        local x
        local selected_models=$(for x in ${available_models} ; do [[ $x == "${cur}"* ]] && echo $x; done)
        COMPREPLY=( $selected_models )
        return 0
    fi

    # handle --with-math-cmd=
    if [[ ${prev} == "--with-math-cmd" && ${cur} == "=" ]] ; then
        COMPREPLY=( $(compgen -c) )
        return 0
    fi

    # handle --with-math-cmd=XXX
    if [[ ${pprev} == "--with-math-cmd" && ${prev} == "=" ]] ; then
        COMPREPLY=( $(compgen -c -- "${cur}") )
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
