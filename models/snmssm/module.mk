DIR          := models/snmssm
MODNAME      := libsnmssm

ifeq ($(shell $(FSCONFIG) --with-smssm),yes)
LIBSNMSSM_SRC  := \
		$(DIR)/nmssmUtils.cpp \
		$(DIR)/nmssmsoftpars.cpp \
		$(DIR)/nmssmsoftsusy.cpp \
		$(DIR)/nmssmsusy.cpp \
		$(DIR)/nmssm1loop.f \
		$(DIR)/nmssm2loop.f

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSNMSSM_SRC  += \
		$(DIR)/snmssm_two_scale.cpp \
		$(DIR)/snmssm_two_scale_convergence_tester.cpp \
		$(DIR)/snmssm_two_scale_initial_guesser.cpp \
		$(DIR)/snmssm_two_scale_low_scale_constraint.cpp \
		$(DIR)/snmssm_two_scale_sugra_constraint.cpp \
		$(DIR)/snmssm_two_scale_susy_scale_constraint.cpp
endif
endif

LIBSNMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSNMSSM_SRC)))

LIBSNMSSM_DEP  := \
		$(LIBSNMSSM_OBJ:.o=.d)

LIBSNMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSNMSSM)

clean-$(MODNAME):
		rm -rf $(LIBSNMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSNMSSM_DEP)
		rm -rf $(LIBSNMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSNMSSM): $(LIBSNMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSNMSSM_DEP)
ALLLIB += $(LIBSNMSSM)
