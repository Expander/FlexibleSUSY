DIR          := models/smssm
MODNAME      := libsmssm

LIBSMSSM_SRC  := \
		$(DIR)/mssmUtils.cpp \
		$(DIR)/physpars.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/twoloophiggs.f

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSMSSM_SRC  += \
		$(DIR)/mssm_two_scale.cpp \
		$(DIR)/mssm_two_scale_convergence_tester.cpp \
		$(DIR)/mssm_two_scale_initial_guesser.cpp \
		$(DIR)/mssm_two_scale_susy_scale_constraint.cpp \
		$(DIR)/mssm_two_scale_low_scale_constraint.cpp \
		$(DIR)/mssm_two_scale_sugra_constraint.cpp
endif

LIBSMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSMSSM_SRC)))

LIBSMSSM_DEP  := \
		$(LIBSMSSM_OBJ:.o=.d)

LIBSMSSM      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSMSSM)

clean-$(MODNAME):
		rm -rf $(LIBSMSSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSMSSM_DEP)
		rm -rf $(LIBSMSSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSMSSM): $(LIBSMSSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSMSSM_DEP)
ALLLIB += $(LIBSMSSM)
