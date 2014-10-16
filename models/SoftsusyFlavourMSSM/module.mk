DIR          := models/SoftsusyFlavourMSSM
MODNAME      := SoftsusyFlavourMSSM

LIBSoftsusyFlavourMSSM_SRC  := \
		$(DIR)/flavoursoft.cpp

LIBSoftsusyFlavourMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyFlavourMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyFlavourMSSM_SRC)))

LIBSoftsusyFlavourMSSM_DEP  := \
		$(LIBSoftsusyFlavourMSSM_OBJ:.o=.d)

LIBSoftsusyFlavourMSSM      := $(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSoftsusyFlavourMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyFlavourMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyFlavourMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyFlavourMSSM)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyFlavourMSSM_DEP) $(LIBSoftsusyFlavourMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBSoftsusyFlavourMSSM): $(LIBSoftsusyFlavourMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSoftsusyFlavourMSSM_DEP)
ALLLIB += $(LIBSoftsusyFlavourMSSM)
