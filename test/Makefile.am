AM_CPPFLAGS = -I$(top_srcdir)/src -I$(BOOST_INC)

bin_PROGRAMS = test_two_scale_solver

test_two_scale_solver_LDADD = \
	$(top_srcdir)/src/libflexisusy.a \
	-lboost_unit_test_framework

test_two_scale_solver_SOURCES = \
	test_two_scale_solver.cpp
