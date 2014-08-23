DIR          := models/sm
MODNAME      := libsm

LIBsm_SRC    :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBsm_SRC    += \
		$(DIR)/sm_two_scale.cpp \
		$(DIR)/sm_two_scale_convergence_tester.cpp \
		$(DIR)/sm_two_scale_experimental_constraint.cpp
endif

LIBsm_OBJ    := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBsm_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBsm_SRC)))

LIBsm_DEP    := \
		$(LIBsm_OBJ:.o=.d)

LIBsm        := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBsm)

clean-$(MODNAME):
		rm -rf $(LIBsm_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBsm_DEP)
		rm -rf $(LIBsm)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBsm): $(LIBsm_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBsm_DEP)
ALLLIB += $(LIBsm)
