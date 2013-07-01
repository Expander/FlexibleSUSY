DIR          := legacy
MODNAME      := liblegacy

LIBLEGACY_SRC := \
		$(DIR)/rk_legacy.cpp

LIBLEGACY_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBLEGACY_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBLEGACY_SRC)))

LIBLEGACY_DEP := \
		$(LIBLEGACY_OBJ:.o=.d)

LIBLEGACY     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBLEGACY)

clean-$(MODNAME):
		rm -rf $(LIBLEGACY_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBLEGACY_DEP)
		rm -rf $(LIBLEGACY)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBLEGACY): $(LIBLEGACY_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBLEGACY_DEP)
ALLLIB += $(LIBLEGACY)
