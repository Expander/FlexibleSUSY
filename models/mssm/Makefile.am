
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(BOOST_INC)

bin_PROGRAMS = run test_MssmSolver

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	susy.cpp \
	softpars.cpp \
	softsusy.cpp \
	physpars.cpp \
	twoloophiggs.f

test_MssmSolver_LDADD = \
	@FLIBS@ \
	$(top_srcdir)/src/libflexisusy.a \
	$(BOOST_LIB)/libboost_unit_test_framework.a

test_MssmSolver_SOURCES = \
	test_MssmSolver.cpp \
	susy.cpp \
	softpars.cpp \
	softsusy.cpp \
	physpars.cpp \
	twoloophiggs.f
