DIR          := higgs-study
MODNAME      := higgs

SCAN_MSSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-MSSM),yes)
SCAN_MSSM_SRC := $(DIR)/scan_MSSM.cpp
endif
SCAN_MSSM_DEP := $(SCAN_MSSM_SRC:.cpp=.d)
SCAN_MSSM_OBJ := $(SCAN_MSSM_SRC:.cpp=.o)
SCAN_MSSM_EXE := $(SCAN_MSSM_SRC:.cpp=.x)

SCAN_MSSM_TC_SRC :=
ifeq ($(shell $(FSCONFIG) --with-MSSMNoFV),yes)
SCAN_MSSM_TC_SRC := $(DIR)/scan_MSSM_tc.cpp
endif
SCAN_MSSM_TC_DEP := $(SCAN_MSSM_TC_SRC:.cpp=.d)
SCAN_MSSM_TC_OBJ := $(SCAN_MSSM_TC_SRC:.cpp=.o)
SCAN_MSSM_TC_EXE := $(SCAN_MSSM_TC_SRC:.cpp=.x)

SCAN_NMSSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-NMSSM),yes)
SCAN_NMSSM_SRC := $(DIR)/scan_NMSSM.cpp
endif
SCAN_NMSSM_DEP := $(SCAN_NMSSM_SRC:.cpp=.d)
SCAN_NMSSM_OBJ := $(SCAN_NMSSM_SRC:.cpp=.o)
SCAN_NMSSM_EXE := $(SCAN_NMSSM_SRC:.cpp=.x)

SCAN_phdUMSSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-phdUMSSM),yes)
SCAN_phdUMSSM_SRC := $(DIR)/scan_phdUMSSM.cpp
endif
SCAN_phdUMSSM_DEP := $(SCAN_phdUMSSM_SRC:.cpp=.d)
SCAN_phdUMSSM_OBJ := $(SCAN_phdUMSSM_SRC:.cpp=.o)
SCAN_phdUMSSM_EXE := $(SCAN_phdUMSSM_SRC:.cpp=.x)

SCAN_phdE6SSM_SRC :=
ifeq ($(shell $(FSCONFIG) --with-phdE6SSM),yes)
SCAN_phdE6SSM_SRC := $(DIR)/scan_phdE6SSM.cpp
endif
SCAN_phdE6SSM_DEP := $(SCAN_phdE6SSM_SRC:.cpp=.d)
SCAN_phdE6SSM_OBJ := $(SCAN_phdE6SSM_SRC:.cpp=.o)
SCAN_phdE6SSM_EXE := $(SCAN_phdE6SSM_SRC:.cpp=.x)

SCAN_NMSSMNoUni_SRC :=
ifeq ($(shell $(FSCONFIG) --with-NMSSMNoUni --with-NMSSM),yes yes)
SCAN_NMSSMNoUni_SRC := $(DIR)/scan_NMSSM_Alambda.cpp
endif
SCAN_NMSSMNoUni_DEP := $(SCAN_NMSSMNoUni_SRC:.cpp=.d)
SCAN_NMSSMNoUni_OBJ := $(SCAN_NMSSMNoUni_SRC:.cpp=.o)
SCAN_NMSSMNoUni_EXE := $(SCAN_NMSSMNoUni_SRC:.cpp=.x)

ALLHIGGS_SRC := \
		$(SCAN_MSSM_SRC) \
		$(SCAN_MSSM_TC_SRC) \
		$(SCAN_NMSSM_SRC) \
		$(SCAN_phdE6SSM_SRC) \
		$(SCAN_phdUMSSM_SRC) \
		$(SCAN_NMSSMNoUni_SRC)

ALLHIGGS_EXE := \
		$(SCAN_MSSM_EXE) \
		$(SCAN_MSSM_TC_EXE) \
		$(SCAN_NMSSM_EXE) \
		$(SCAN_phdE6SSM_EXE) \
		$(SCAN_phdUMSSM_EXE) \
		$(SCAN_NMSSMNoUni_EXE)

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

$(SCAN_MSSM_TC_EXE): $(SCAN_MSSM_TC_OBJ) $(LIBMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(SCAN_NMSSM_EXE): $(SCAN_NMSSM_OBJ) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(SCAN_phdE6SSM_EXE): $(SCAN_phdE6SSM_OBJ) $(LIBphdE6SSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(SCAN_phdUMSSM_EXE): $(SCAN_phdUMSSM_OBJ) $(LIBphdUMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(SCAN_NMSSMNoUni_EXE): $(SCAN_NMSSMNoUni_OBJ) $(LIBNMSSM) $(LIBNMSSMNoUni) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

ALLDEP += $(ALLHIGGS_DEP)
ALLEXE += $(ALLHIGGS_EXE)
