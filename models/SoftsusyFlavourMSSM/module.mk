DIR          := models/SoftsusyFlavourMSSM
MODNAME      := SoftsusyFlavourMSSM
WITH_$(MODNAME) := yes
MODSoftsusyFlavourMSSM_MOD := MSSM_higgs NMSSM_higgs
MODSoftsusyFlavourMSSM_DEP := $(patsubst %,model_specific/%,$(MODSoftsusyFlavourMSSM_MOD))
MODSoftsusyFlavourMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSoftsusyFlavourMSSM_MOD))
MODSoftsusyFlavourMSSM_LIB := $(foreach M,$(MODSoftsusyFlavourMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

LIBSoftsusyFlavourMSSM_SRC  := \

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

$(LIBSoftsusyFlavourMSSM_DEP) $(LIBSoftsusyFlavourMSSM_OBJ) $(EXESoftsusyFlavourMSSM_DEP) $(EXESoftsusyFlavourMSSM_OBJ): \
	CPPFLAGS += $(MODSoftsusyFlavourMSSM_INC) $(EIGENFLAGS)

$(LIBSoftsusyFlavourMSSM): $(LIBSoftsusyFlavourMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBSoftsusyFlavourMSSM_DEP) $(EXESoftsusyFlavourMSSM_DEP)
ALLLIB += $(LIBSoftsusyFlavourMSSM)
ALLEXE += $(EXESoftsusyFlavourMSSM_EXE)
ALLMODDEP += $(MODSoftsusyFlavourMSSM_DEP)
