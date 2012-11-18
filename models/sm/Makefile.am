AM_CPPFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/models/sm \
	@BOOST_CPPFLAGS@

bin_PROGRAMS = run

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	sm_two_scale.cpp \
	sm_two_scale_experimental_constraint.cpp
