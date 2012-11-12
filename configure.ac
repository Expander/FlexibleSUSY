AC_INIT(SUSY Spectrum Generator,0.1,Alexander.Voigt@physik.tu-dresden.de,FlexibleSUSY)
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR(src/Makefile.am)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)

sysconfdir=${sysconfdir}/${PACKAGE}/${VERSION}

AC_PROG_CXX
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_PROG_RANLIB

m4_include([m4/m4_ax_boost_base.m4])
m4_include([m4/m4_ax_boost_unit_test_framework.m4])
AX_BOOST_BASE([1.33.1])
AX_BOOST_UNIT_TEST_FRAMEWORK

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/test/unit_test.hpp], [],
    [AC_MSG_ERROR(You need the Boost libraries.)])
AC_LANG_POP([C++])

AM_CXXFLAGS="-ggdb -Wall -pedantic -Wextra -Wcast-qual -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor -Wfloat-equal"

AC_SUBST(FLIBS)
AC_SUBST(AM_CXXFLAGS)

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([test/Makefile])
AC_CONFIG_FILES([models/Makefile])
AC_CONFIG_FILES([models/sm/Makefile])
AC_CONFIG_FILES([models/mssm/Makefile])

AC_OUTPUT
