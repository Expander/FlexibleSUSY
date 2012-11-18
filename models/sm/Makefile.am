AM_CPPFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/models/sm \
	@BOOST_CPPFLAGS@

bin_PROGRAMS = \
	run \
	test_sm_two_scale

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	sm_two_scale.cpp \
	sm_two_scale_experimental_constraint.cpp

test_sm_two_scale_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	@BOOST_LDFLAGS@ \
	@BOOST_UNIT_TEST_FRAMEWORK_LIB@

test_sm_two_scale_SOURCES = \
	test_sm_two_scale.cpp \
	sm_two_scale.cpp \
	sm_two_scale_experimental_constraint.cpp
