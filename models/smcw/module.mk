DIR          := models/smcw
MODNAME      := libsmcw

LIBSMCW_HDR  := \
		$(DIR)/smcw.hpp \
		$(DIR)/smcw_two_scale.hpp \
		$(DIR)/smcw_two_scale_convergence_tester.hpp \
		$(DIR)/smcw_two_scale_gut_constraint.hpp

LIBSMCW_SRC  := \
		$(DIR)/smcw_two_scale.cpp \
		$(DIR)/smcw_two_scale_convergence_tester.cpp \
		$(DIR)/smcw_two_scale_gut_constraint.cpp

LIBSMCW_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSMCW_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSMCW_SRC)))

LIBSMCW_DEP  := \
		$(LIBSMCW_OBJ:.o=.d)

LIBSMCW      := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSMCW)

clean-$(MODNAME):
		rm -rf $(LIBSMCW_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSMCW_DEP)
		rm -rf $(LIBSMCW)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSMCW): $(LIBSMCW_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSMCW_DEP)
ALLLIB += $(LIBSMCW)
