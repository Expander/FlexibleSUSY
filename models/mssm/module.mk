DIR          := models/mssm
MODNAME      := libmssm

LIBMSSM_SRC  := \
		$(DIR)/physpars.cpp \
		$(DIR)/softpars.cpp \
		$(DIR)/softsusy.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/twoloophiggs.f

ifneq ($(findstring two_scale,$(ALGORITMS)),)
LIBMSSM_SRC  += \
		$(DIR)/mssm_two_scale.cpp \
		$(DIR)/mssm_two_scale_convergence_tester.cpp \
		$(DIR)/mssm_two_scale_initial_guesser.cpp \
		$(DIR)/mssm_two_scale_msusy_constraint.cpp \
		$(DIR)/mssm_two_scale_mz_constraint.cpp \
		$(DIR)/mssm_two_scale_sugra_constraint.cpp
endif

LIBMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSM_SRC)))

LIBMSSM_DEP  := \
		$(LIBMSSM_OBJ:.o=.d)

LIBMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBMSSM)

clean-$(MODNAME):
		rm -rf $(LIBMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBMSSM_DEP)
		rm -rf $(LIBMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBMSSM): $(LIBMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBMSSM_DEP)
ALLLIB += $(LIBMSSM)
