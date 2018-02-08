DIR          := models/SoftsusyFlavourMSSM
MODNAME      := SoftsusyFlavourMSSM
WITH_$(MODNAME) := yes

LIBSoftsusyFlavourMSSM_SRC  := \
		$(DIR)/flavoursoft.cpp

LIBSoftsusyFlavourMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyFlavourMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyFlavourMSSM_SRC)))

LIBSoftsusyFlavourMSSM_DEP  := \
		$(LIBSoftsusyFlavourMSSM_OBJ:.o=.d)

LIBSoftsusyFlavourMSSM      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

EXESoftsusyFlavourMSSM_SRC  := $(DIR)/run_softpoint.cpp

EXESoftsusyFlavourMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESoftsusyFlavourMSSM_SRC)))

EXESoftsusyFlavourMSSM_DEP  := \
		$(EXESoftsusyFlavourMSSM_OBJ:.o=.d)

EXESoftsusyFlavourMSSM_EXE := \
		$(EXESoftsusyFlavourMSSM_OBJ:.o=.x)

RUN_SOFTPOINT_EXE := $(EXESoftsusyFlavourMSSM_EXE)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSoftsusyFlavourMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyFlavourMSSM_DEP)
		-rm -f $(EXESoftsusyFlavourMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyFlavourMSSM_OBJ)
		-rm -f $(EXESoftsusyFlavourMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyFlavourMSSM)
		-rm -f $(EXESoftsusyFlavourMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyFlavourMSSM_DEP) $(LIBSoftsusyFlavourMSSM_OBJ) $(EXESoftsusyFlavourMSSM_DEP) $(EXESoftsusyFlavourMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBSoftsusyFlavourMSSM): $(LIBSoftsusyFlavourMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(EXESoftsusyFlavourMSSM_EXE): $(DIR)/run_softpoint.o $(LIBSoftsusyFlavourMSSM) $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIB_model_specific_MSSM_higgs) $(LIB_model_specific_NMSSM_higgs) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

ALLDEP += $(LIBSoftsusyFlavourMSSM_DEP) $(EXESoftsusyFlavourMSSM_DEP)
ALLLIB += $(LIBSoftsusyFlavourMSSM)
ALLEXE += $(EXESoftsusyFlavourMSSM_EXE)
