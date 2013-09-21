DIR          := models/SoftsusyNMSSM
MODNAME      := libSoftsusyNMSSM

ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM),yes)
LIBSoftsusyNMSSM_SRC  := \
		$(DIR)/nmssmUtils.cpp \
		$(DIR)/nmssmsoftpars.cpp \
		$(DIR)/nmssmsoftsusy.cpp \
		$(DIR)/nmssmsusy.cpp \
		$(DIR)/nmssm1loop.f \
		$(DIR)/nmssm2loop.f

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSoftsusyNMSSM_SRC  += \
		$(DIR)/snmssm_two_scale.cpp \
		$(DIR)/snmssm_two_scale_convergence_tester.cpp \
		$(DIR)/snmssm_two_scale_initial_guesser.cpp \
		$(DIR)/snmssm_two_scale_low_scale_constraint.cpp \
		$(DIR)/snmssm_two_scale_sugra_constraint.cpp \
		$(DIR)/snmssm_two_scale_susy_scale_constraint.cpp
endif
endif

LIBSoftsusyNMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyNMSSM_SRC)))

LIBSoftsusyNMSSM_DEP  := \
		$(LIBSoftsusyNMSSM_OBJ:.o=.d)

LIBSoftsusyNMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSoftsusyNMSSM)

clean-$(MODNAME):
		rm -rf $(LIBSoftsusyNMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSoftsusyNMSSM_DEP)
		rm -rf $(LIBSoftsusyNMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyNMSSM): $(LIBSoftsusyNMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSoftsusyNMSSM_DEP)
ALLLIB += $(LIBSoftsusyNMSSM)
