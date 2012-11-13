AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/models/sm

bin_PROGRAMS = run

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	smcw_two_scale.cpp \
	$(top_srcdir)/models/sm/sm_two_scale.cpp
