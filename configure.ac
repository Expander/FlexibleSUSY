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

AC_SUBST(FLIBS)
AC_SUBST([AM_CXXFLAGS], [-ggdb -Wall -pedantic -Wextra -Wcast-qual -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor -Wfloat-equal])

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([test/Makefile])
AC_CONFIG_FILES([models/Makefile])
AC_CONFIG_FILES([models/sm/Makefile])

AC_OUTPUT
