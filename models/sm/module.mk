DIR          := models/sm
MODNAME      := libsm

LIBSM_SRC    :=

ifneq ($(findstring two_scale,$(ALGORITMS)),)
LIBSM_SRC    += \
		$(DIR)/sm_two_scale.cpp \
		$(DIR)/sm_two_scale_convergence_tester.cpp \
		$(DIR)/sm_two_scale_experimental_constraint.cpp
endif

LIBSM_OBJ    := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSM_SRC)))

LIBSM_DEP    := \
		$(LIBSM_OBJ:.o=.d)

LIBSM        := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSM)

clean-$(MODNAME):
		rm -rf $(LIBSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSM_DEP)
		rm -rf $(LIBSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSM): $(LIBSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSM_DEP)
ALLLIB += $(LIBSM)
