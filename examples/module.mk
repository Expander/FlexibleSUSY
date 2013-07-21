DIR      := examples
MODNAME  := examples

ifeq ($(shell $(FSCONFIG) --with-smssm),yes)
EXAMPLES_SRC := \
		$(DIR)/run_softsusy.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
EXAMPLES_SRC += \
		$(DIR)/run_smssm.cpp
endif
endif

ifneq ($(findstring lattice,$(ALGORITHMS)),)
ifeq ($(shell $(FSCONFIG) --with-fmssm),yes)
LATTICE_EXAMPLES_SRC += \
		$(DIR)/lattice_fmssm.cpp \
		$(DIR)/lattice_numerical_fmssm.cpp \
		$(DIR)/lattice_fmssm_fmssmn.cpp \
		$(DIR)/lattice_numerical_fmssm_fmssmn.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-fmssm --with-fmssn),yes yes)
LATTICE_EXAMPLES_SRC += \
		$(DIR)/lattice_fmssm_fmssmn.cpp \
		$(DIR)/lattice_numerical_fmssm_fmssmn.cpp
endif

LATTICE_EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LATTICE_EXAMPLES_SRC)))

LATTICE_EXAMPLES_DEP := \
		$(LATTICE_EXAMPLES_OBJ:.o=.d)

EXAMPLES_SRC += $(LATTICE_EXAMPLES_SRC)
endif

EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXAMPLES_SRC)))

EXAMPLES_DEP := \
		$(EXAMPLES_OBJ:.o=.d)

EXAMPLES_EXE := \
		$(EXAMPLES_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(EXAMPLES_EXE)

clean-$(MODNAME):
		rm -rf $(EXAMPLES_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(EXAMPLES_DEP)
		rm -rf $(EXAMPLES_EXE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/run_smssm.d $(DIR)/run_smssm.o: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/run_smssm.x: $(DIR)/run_smssm.o $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(FLIBS)

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LATTICE_EXAMPLES_DEP) $(LATTICE_EXAMPLES_OBJ): CPPFLAGS += $(EIGENFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)

$(DIR)/lattice_fmssm.x: $(DIR)/lattice_fmssm.o $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS)

$(DIR)/lattice_numerical_fmssm.x: $(DIR)/lattice_numerical_fmssm.o $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS)

$(DIR)/lattice_fmssm_fmssmn.x: $(DIR)/lattice_fmssm_fmssmn.o \
			       $(LIBFMSSMN) $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS)
$(DIR)/lattice_numerical_fmssm_fmssmn.x: $(DIR)/lattice_numerical_fmssm_fmssmn.o \
			       $(LIBFMSSMN) $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS)
endif

$(DIR)/run_softsusy.x: $(DIR)/run_softsusy.o $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(abspath $^) $(FLIBS)

ALLDEP += $(EXAMPLES_DEP)
ALLEXE += $(EXAMPLES_EXE)
