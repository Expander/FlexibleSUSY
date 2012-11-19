MODULES  := src test
CXXFLAGS := $(patsubst %, -I%, $(MODULES))
FFLAGS   := $(patsubst %, -I%, $(MODULES))
LIBS     := -lboost_unit_test_framework
ALLDEP   :=
ALLLIB   :=
ALLEXE   :=
CXX      := g++
FC       := gfortran
MAKELIB  := ar r

.PHONY:  all clean distclean

all:     allexec alllib

include $(patsubst %, %/module.mk, $(MODULES))

ifeq ($(findstring $(MAKECMDGOALS),clean distclean),)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
-include $(ALLDEP)
endif
endif

allexec:  $(ALLEXE)
alllib:   $(ALLLIB)

%.d: %.cpp
# -MT '$*.o' ensures that the target contains the full path
	$(CXX) $(CXXFLAGS) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f
# the sed script ensures that the target contains the full path
	$(FC) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

clean::
	rm -f $(ALLDEP)

distclean:: clean
	rm -f $(ALLLIB) $(ALLEXE)
