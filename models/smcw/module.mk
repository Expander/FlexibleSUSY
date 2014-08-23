DIR          := models/smcw
MODNAME      := libsmcw

LIBsmcw_SRC  :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBsmcw_SRC  += \
		$(DIR)/smcw_two_scale.cpp \
		$(DIR)/smcw_two_scale_convergence_tester.cpp \
		$(DIR)/smcw_two_scale_gut_constraint.cpp
endif

LIBsmcw_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBsmcw_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBsmcw_SRC)))

LIBsmcw_DEP  := \
		$(LIBsmcw_OBJ:.o=.d)

LIBsmcw      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBsmcw)

clean-$(MODNAME):
		rm -rf $(LIBsmcw_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBsmcw_DEP)
		rm -rf $(LIBsmcw)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBsmcw): $(LIBsmcw_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBsmcw_DEP)
ALLLIB += $(LIBsmcw)
