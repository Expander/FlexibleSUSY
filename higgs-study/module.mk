DIR          := higgs-study
MODNAME      := higgs

SCAN_MSSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-MSSM),yes)
SCAN_MSSM_SRC := $(DIR)/scanMSSM.cpp
endif
SCAN_MSSM_DEP := $(SCAN_MSSM_SRC:.cpp=.d)
SCAN_MSSM_OBJ := $(SCAN_MSSM_SRC:.cpp=.o)
SCAN_MSSM_EXE := $(SCAN_MSSM_SRC:.cpp=.x)

SCAN_NMSSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-NMSSM),yes)
SCAN_NMSSM_SRC := $(DIR)/scanNMSSM.cpp
endif
SCAN_NMSSM_DEP := $(SCAN_NMSSM_SRC:.cpp=.d)
SCAN_NMSSM_OBJ := $(SCAN_NMSSM_SRC:.cpp=.o)
SCAN_NMSSM_EXE := $(SCAN_NMSSM_SRC:.cpp=.x)

ALLHIGGS_SRC := \
		$(SCAN_MSSM_SRC) \
		$(SCAN_NMSSM_SRC)

ALLHIGGS_EXE := \
		$(SCAN_MSSM_EXE) \
		$(SCAN_NMSSM_EXE)

ALLHIGGS_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(ALLHIGGS_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(ALLHIGGS_SRC)))

ALLHIGGS_DEP := \
		$(ALLHIGGS_OBJ:.o=.d)

ALLHIGGS     := $(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(ALLHIGGS_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(ALLHIGGS_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(ALLHIGGS_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(ALLHIGGS)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(ALLHIGGS_DEP) $(ALLHIGGS_OBJ): CPPFLAGS += $(EIGENFLAGS)
$(ALLHIGGS_DEP) $(ALLHIGGS_OBJ): CPPFLAGS += $(GSLFLAGS) $(BOOSTFLAGS)

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(ALLHIGGS): $(ALLHIGGS_OBJ)
		$(MAKELIB) $@ $^
else
$(ALLHIGGS): $(ALLHIGGS_OBJ)
		$(MAKELIB) $@ $^ $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS)
endif

$(SCAN_MSSM_EXE): $(SCAN_MSSM_OBJ) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(SCAN_NMSSM_EXE): $(SCAN_NMSSM_OBJ) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

ALLDEP += $(ALLHIGGS_DEP)
ALLEXE += $(ALLHIGGS_EXE)
