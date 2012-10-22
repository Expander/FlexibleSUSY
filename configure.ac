AC_INIT(SUSY Spectrum Generator,0.1,Alexander.Voigt@physik.tu-dresden.de,FlexibleSUSY)
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR(src/main.cpp)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)

sysconfdir=${sysconfdir}/${PACKAGE}/${VERSION}

AC_PROG_CXX
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

AC_SUBST(FLIBS)

AC_OUTPUT(Makefile src/Makefile doc/Makefile test/Makefile)
