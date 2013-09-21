DIR          := models/SoftsusyMSSM
MODNAME      := libSoftsusyMSSM

LIBSoftsusyMSSM_SRC  := \
		$(DIR)/mssmUtils.cpp \
		$(DIR)/physpars.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/twoloophiggs.f

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSoftsusyMSSM_SRC  += \
		$(DIR)/mssm_two_scale.cpp \
		$(DIR)/mssm_two_scale_convergence_tester.cpp \
		$(DIR)/mssm_two_scale_initial_guesser.cpp \
		$(DIR)/mssm_two_scale_susy_scale_constraint.cpp \
		$(DIR)/mssm_two_scale_low_scale_constraint.cpp \
		$(DIR)/mssm_two_scale_sugra_constraint.cpp
endif

LIBSoftsusyMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyMSSM_SRC)))

LIBSoftsusyMSSM_DEP  := \
		$(LIBSoftsusyMSSM_OBJ:.o=.d)

LIBSoftsusyMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSoftsusyMSSM)

clean-$(MODNAME):
		rm -rf $(LIBSoftsusyMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSoftsusyMSSM_DEP)
		rm -rf $(LIBSoftsusyMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyMSSM): $(LIBSoftsusyMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSoftsusyMSSM_DEP)
ALLLIB += $(LIBSoftsusyMSSM)
