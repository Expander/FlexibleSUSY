DIR          := models/SoftsusyNMSSM
MODNAME      := SoftsusyNMSSM

ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM),yes)
LIBSoftsusyNMSSM_SRC  := \
		$(DIR)/nmssmUtils.cpp \
		$(DIR)/nmssmsoftpars.cpp \
		$(DIR)/nmssmsoftsusy.cpp \
		$(DIR)/nmssmsusy.cpp \
		$(DIR)/nmssm1loop.f

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSoftsusyNMSSM_SRC  += \
		$(DIR)/SoftsusyNMSSM_two_scale.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_sugra_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_susy_scale_constraint.cpp
endif
endif

LIBSoftsusyNMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyNMSSM_SRC)))

LIBSoftsusyNMSSM_DEP  := \
		$(LIBSoftsusyNMSSM_OBJ:.o=.d)

LIBSoftsusyNMSSM      := $(DIR)/lib$(MODNAME)$(LIBEXT)

EXESoftsusyNMSSM_SRC  :=

ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-SoftsusyNMSSM),yes yes)
EXESoftsusyNMSSM_SRC  += \
		$(DIR)/run_softpoint.cpp
endif

EXESoftsusyNMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESoftsusyNMSSM_SRC)))

EXESoftsusyNMSSM_DEP  := \
		$(EXESoftsusyNMSSM_OBJ:.o=.d)

RUN_SOFTPOINT_EXE := \
		$(EXESoftsusyNMSSM_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj

all-$(MODNAME): $(LIBSoftsusyNMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyNMSSM_DEP)
		-rm -f $(EXESoftsusyNMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyNMSSM_OBJ)
		-rm -f $(EXESoftsusyNMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyNMSSM)
		-rm -f $(RUN_SOFTPOINT_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyNMSSM): $(LIBSoftsusyNMSSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_SOFTPOINT_EXE): $(DIR)/run_softpoint.o $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

ALLDEP += $(LIBSoftsusyNMSSM_DEP) $(EXESoftsusyNMSSM_DEP)
ALLLIB += $(LIBSoftsusyNMSSM)
ALLEXE += $(RUN_SOFTPOINT_EXE)
