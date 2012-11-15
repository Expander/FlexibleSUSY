AM_CPPFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/models/sm \
	-I$(top_srcdir)/models/smcw \
	@BOOST_CPPFLAGS@

bin_PROGRAMS = \
	test_two_scale_solver \
	test_sm_smcw_matching \
	test_logger

test_two_scale_solver_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	@BOOST_LDFLAGS@ \
	@BOOST_UNIT_TEST_FRAMEWORK_LIB@

test_two_scale_solver_SOURCES = \
	test_two_scale_solver.cpp

test_sm_smcw_matching_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	@BOOST_LDFLAGS@ \
	@BOOST_UNIT_TEST_FRAMEWORK_LIB@

test_sm_smcw_matching_SOURCES = \
	test_sm_smcw_matching.cpp \
	$(top_srcdir)/models/sm/sm_two_scale.cpp \
	$(top_srcdir)/models/smcw/smcw_two_scale.cpp

test_logger_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	@BOOST_LDFLAGS@ \
	@BOOST_UNIT_TEST_FRAMEWORK_LIB@

test_logger_SOURCES = \
	test_logger.cpp
