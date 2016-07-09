#! /bin/sh

# This script returns the include paths for Mathematica's MathLink and
# LibraryLink .
# Author: Alexander Voigt

MATH_CMD=${MATH:-math}

eval `"${MATH_CMD}" -run '
    Print["sysid=\"", $SystemID, "\""];
    Print["topdir=\"", $TopDirectory, "\""];
    Exit[]' < /dev/null | tr '\r' ' ' | tail -2`

# check whether Cygwin's dlltool can handle 64-bit DLLs
test "$sysid" = Windows-x86-64 && {
    ${DLLTOOL:-dlltool} --help | grep x86-64 > /dev/null || sysid=Windows
}

topdir=`cd "$topdir" ; echo $PWD`

get_librarylink_incpath() {
    for p in \
        "$topdir/SystemFiles/IncludeFiles/C" ; do
        test -d "$p" && break
    done

    echo "$p"
}

get_mathlink_incpath() {
    for p in \
        "$topdir/SystemFiles/Links/MathLink/DeveloperKit/$sysid/CompilerAdditions" \
        "$topdir/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions" \
        "$topdir/AddOns/MathLink/DeveloperKit/$sysid/CompilerAdditions" ; do
        test -d "$p" && break
    done

    echo "$p"
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case $1 in
            -I)               prepend="-I" ;;
            --librarylink|-l) path="${path} '${prepend}$(get_librarylink_incpath)'" ;;
            --mathlink|-m)    path="${path} '${prepend}$(get_mathlink_incpath)'" ;;
            --help|-h)        usage; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

echo "$path"
