MODULES  := config src models/sm models/smcw models/mssm test doc

BOOSTLIBDIR :=
BOOSTLIBS   := -lboost_unit_test_framework
BOOSTINC    :=
ifneq "$(BOOSTLIBDIR)" ""
   BOOSTLIBS := $(addprefix -L$(BOOSTLIBDIR) , $(BOOSTLIBS))
endif

CPPFLAGS := $(BOOSTINC) $(patsubst %, -I%, $(MODULES))
CXXFLAGS := -ggdb -Wall -pedantic -Wextra -Wcast-qual \
            -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor \
            -g -O2
FFLAGS   :=
LIBS     := $(BOOSTLIBS) -lgfortranbegin -lgfortran -lm
CXX      := g++
FC       := gfortran
MAKELIB  := ar cru
VERSION  := 0.1
PKGNAME  := FlexibleSUSY

# the modules add their dependency files to this variable
ALLDEP   :=
# the modules add headers to be created to this variable
ALLHDR   :=
# the modules add their libraries to this variable
ALLLIB   :=
# the modules add executables to this variable
ALLEXE   :=

.PHONY:  all allhdr allexec alllib clean distclean tag release

all:     allhdr allexec alllib

include $(patsubst %, %/module.mk, $(MODULES))

ifeq ($(findstring $(MAKECMDGOALS),clean distclean tag release),)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
-include $(ALLDEP)
endif
endif
endif

allhdr:   $(ALLHDR)
allexec:  $(ALLEXE)
alllib:   $(ALLLIB)

%.d: %.cpp
# -MT '$*.o' ensures that the target contains the full path
	$(CXX) $(CPPFLAGS) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f
# the sed script ensures that the target contains the full path
	$(FC) $(CPPFLAGS) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

tag:
	git tag v$(VERSION) -m "version $(VERSION)"

release:
	git archive --worktree-attributes --format=tar \
	--prefix=$(PKGNAME)-$(VERSION)/ \
	v$(VERSION) | gzip > $(PKGNAME)-$(VERSION).tar.gz
