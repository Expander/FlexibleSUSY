DIR      := examples
MODNAME  := examples

EXAMPLES_SRC :=

ifneq ($(findstring lattice,$(ALGORITHMS)),)
ifeq ($(shell $(FSCONFIG) --with-fmssm),yes)
LATTICE_EXAMPLES_SRC := \
		$(DIR)/lattice_fmssm.cpp \
		$(DIR)/lattice_numerical_fmssm.cpp
ifeq ($(shell $(FSCONFIG) --with-fmssmn),yes)
LATTICE_EXAMPLES_SRC += \
		$(DIR)/lattice_fmssm_fmssmn.cpp \
		$(DIR)/lattice_numerical_fmssm_fmssmn.cpp
endif
endif

LATTICE_EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LATTICE_EXAMPLES_SRC)))

LATTICE_EXAMPLES_DEP := \
		$(LATTICE_EXAMPLES_OBJ:.o=.d)

EXAMPLES_SRC += $(LATTICE_EXAMPLES_SRC)
endif

ifdef BUILD_SWITCH_EXAMPLES
ifneq ($(findstring two_scale,$(ALGORITHMS)),)
ifneq ($(findstring lattice,$(ALGORITHMS)),)
ifeq ($(shell $(FSCONFIG) --with-fmssm --with-MSSM),yes yes)
SWITCH_EXAMPLES_SRC := \
		$(DIR)/switch_MSSM.cpp

SWITCH_EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(SWITCH_EXAMPLES_SRC)))

SWITCH_EXAMPLES_DEP := \
		$(SWITCH_EXAMPLES_OBJ:.o=.d)

EXAMPLES_SRC += $(SWITCH_EXAMPLES_SRC)
endif
endif
endif
endif

EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXAMPLES_SRC)))

EXAMPLES_DEP := \
		$(EXAMPLES_OBJ:.o=.d)

EXAMPLES_EXE := \
		$(EXAMPLES_OBJ:.o=.x)

STANDALONE_DIR := \
		$(DIR)/standalone-model \
		$(DIR)/standalone-rge

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(EXAMPLES_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(EXAMPLES_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(EXAMPLES_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(EXAMPLES_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		-@for d in $(STANDALONE_DIR); do \
			(cd $$d && make distclean); \
		 done

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LATTICE_EXAMPLES_DEP) $(LATTICE_EXAMPLES_OBJ): CPPFLAGS += $(EIGENFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)

$(DIR)/lattice_fmssm.x: $(DIR)/lattice_fmssm.o $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)

$(DIR)/lattice_numerical_fmssm.x: $(DIR)/lattice_numerical_fmssm.o $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)

$(DIR)/lattice_fmssm_fmssmn.x: $(DIR)/lattice_fmssm_fmssmn.o \
			       $(LIBFMSSMN) $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)
$(DIR)/lattice_numerical_fmssm_fmssmn.x: $(DIR)/lattice_numerical_fmssm_fmssmn.o \
			       $(LIBFMSSMN) $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)
endif

ifdef BUILD_SWITCH_EXAMPLES
ifneq ($(findstring two_scale,$(ALGORITHMS)),)
ifneq ($(findstring lattice,$(ALGORITHMS)),)
ifeq ($(shell $(FSCONFIG) --with-fmssm --with-MSSM),yes yes)
$(SWITCH_EXAMPLES_DEP) $(SWITCH_EXAMPLES_OBJ): CPPFLAGS += $(EIGENFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)

$(DIR)/switch_MSSM.x: $(DIR)/switch_MSSM.o $(LIBMSSM) $(LIBFMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $(call abspathx,$^) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(LOOPTOOLSLIBS) $(LIBFFLITE) $(FLIBS)
endif
endif
endif
endif

ALLDEP += $(EXAMPLES_DEP)
ALLEXE += $(EXAMPLES_EXE)
