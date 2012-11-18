AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/models/sm

EXTRA_DIST = \
	smcw.hpp \
	smcw_two_scale.hpp

bin_PROGRAMS = run

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	smcw_two_scale.cpp \
	$(top_srcdir)/models/sm/sm_two_scale.cpp
