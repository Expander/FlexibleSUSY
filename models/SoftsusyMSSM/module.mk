DIR          := models/SoftsusyMSSM
MODNAME      := SoftsusyMSSM

LIBSoftsusyMSSM_SRC  := \
		$(DIR)/mssmUtils.cpp \
		$(DIR)/physpars.cpp \
		$(DIR)/susy.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSoftsusyMSSM_SRC  += \
		$(DIR)/SoftsusyMSSM_two_scale.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_susy_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_sugra_constraint.cpp
endif

LIBSoftsusyMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyMSSM_SRC)))

LIBSoftsusyMSSM_DEP  := \
		$(LIBSoftsusyMSSM_OBJ:.o=.d)

LIBSoftsusyMSSM      := $(DIR)/lib$(MODNAME)$(LIBEXT)

EXESoftsusyMSSM_SRC := \
		$(DIR)/run_SoftsusyMSSM.cpp

EXESoftsusyMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESoftsusyMSSM_SRC)))

EXESoftsusyMSSM_DEP := \
		$(EXESoftsusyMSSM_OBJ:.o=.d)

RUN_SoftsusyMSSM_EXE := \
		$(EXESoftsusyMSSM_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSoftsusyMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyMSSM_DEP)
		-rm -f $(EXESoftsusyMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyMSSM_OBJ)
		-rm -f $(EXESoftsusyMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyMSSM)
		-rm -f $(RUN_SoftsusyMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyMSSM_DEP) $(EXESoftsusyMSSM_DEP) $(LIBSoftsusyMSSM_OBJ) $(EXESoftsusyMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBSoftsusyMSSM): $(LIBSoftsusyMSSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_SoftsusyMSSM_EXE): $(EXESoftsusyMSSM_OBJ) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

ALLDEP += $(LIBSoftsusyMSSM_DEP) $(EXESoftsusyMSSM_DEP)
ALLLIB += $(LIBSoftsusyMSSM)
ALLEXE += $(RUN_SoftsusyMSSM_EXE)
