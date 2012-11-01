SUBDIRS = test

AM_CPPFLAGS = -I$(top_srcdir)/src

bin_PROGRAMS = run

run_LDADD = @FLIBS@ $(top_srcdir)/src/libflexisusy.a

run_SOURCES = \
	main.cpp \
	susy.cpp \
	softpars.cpp \
	softsusy.cpp \
	physpars.cpp \
	twoloophiggs.f
