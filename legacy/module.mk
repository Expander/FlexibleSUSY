DIR          := legacy
MODNAME      := liblegacy

LIBLEGACY_SRC := \
		$(DIR)/conversion.cpp \
		$(DIR)/diagonalization.cpp \
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
		-rm -f $(LIBLEGACY_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(LIBLEGACY_DEP)
		-rm -f $(LIBLEGACY)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBLEGACY): $(LIBLEGACY_OBJ)
		$(MAKELIB) $@ $^

# add boost and eigen flags for the test object files and dependencies
$(LIBLEGACY_OBJ) $(LIBLEGACY_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ALLDEP += $(LIBLEGACY_DEP)
ALLLIB += $(LIBLEGACY)
