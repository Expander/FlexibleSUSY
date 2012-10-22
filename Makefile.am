SUBDIRS = src doc test
EXTRA_DIST = autogen.sh

datadir = $(datarootdir)/@PACKAGE@/@VERSION@

# create configmake.h, which contains absolute platform install paths
configmake.h: Makefile
	rm -f $@-t $@
	{ echo '/* DO NOT EDIT! GENERATED AUTOMATICALLY! */'; \
	  echo '#define PREFIX "$(prefix)/"'; \
	  echo '#define EXEC_PREFIX "$(exec_prefix)/"'; \
	  echo '#define BINDIR "$(bindir)/"'; \
	  echo '#define SBINDIR "$(sbindir)/"'; \
	  echo '#define LIBEXECDIR "$(libexecdir)/"'; \
	  echo '#define DATAROOTDIR "$(datarootdir)/"'; \
	  echo '#define DATADIR "$(datadir)/"'; \
	  echo '#define SYSCONFDIR "$(sysconfdir)/"'; \
	  echo '#define SHAREDSTATEDIR "$(sharedstatedir)/"'; \
	  echo '#define LOCALSTATEDIR "$(localstatedir)/"'; \
	  echo '#define INCLUDEDIR "$(includedir)/"'; \
	  echo '#define OLDINCLUDEDIR "$(oldincludedir)/"'; \
	  echo '#define DOCDIR "$(docdir)/"'; \
	  echo '#define INFODIR "$(infodir)/"'; \
	  echo '#define HTMLDIR "$(htmldir)/"'; \
	  echo '#define DVIDIR "$(dvidir)/"'; \
	  echo '#define PDFDIR "$(pdfdir)/"'; \
	  echo '#define PSDIR "$(psdir)/"'; \
	  echo '#define LIBDIR "$(libdir)/"'; \
	  echo '#define LISPDIR "$(lispdir)/"'; \
	  echo '#define LOCALEDIR "$(localedir)/"'; \
	  echo '#define MANDIR "$(mandir)/"'; \
	  echo '#define MANEXT "$(manext)/"'; \
	  echo '#define PKGDATADIR "$(pkgdatadir)/"'; \
	  echo '#define PKGINCLUDEDIR "$(pkgincludedir)/"'; \
	  echo '#define PKGLIBDIR "$(pkglibdir)/"'; \
	  echo '#define PKGLIBEXECDIR "$(pkglibexecdir)/"'; \
	} | sed '/""/d' > $@-t
	mv $@-t $@

BUILT_SOURCES = configmake.h
CLEANFILES = configmake.h configmake.h-t
