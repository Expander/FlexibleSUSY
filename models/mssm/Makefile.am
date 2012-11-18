
AM_CPPFLAGS = -I$(top_srcdir)/src

EXTRA_DIST = \
	mssm_solver.h \
	softpars.h \
	susy.h \
	physpars.h \
	softsusy.h \
	twoloophiggs.h

bin_PROGRAMS = run

run_LDADD = \
	@FLIBS@ \
	$(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	susy.cpp \
	softpars.cpp \
	softsusy.cpp \
	physpars.cpp \
	twoloophiggs.f
