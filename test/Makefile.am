AM_CPPFLAGS = -I$(top_srcdir)/src @BOOST_CPPFLAGS@

bin_PROGRAMS = test_two_scale_solver

test_two_scale_solver_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	@BOOST_LDFLAGS@ \
	@BOOST_UNIT_TEST_FRAMEWORK_LIB@

test_two_scale_solver_SOURCES = \
	test_two_scale_solver.cpp
