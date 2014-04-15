DIR          := legacy
MODNAME      := legacy

LIBLEGACY_SRC := \
		$(DIR)/conversion.cpp \
		$(DIR)/diagonalization.cpp \
		$(DIR)/rk_legacy.cpp

LIBLEGACY_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBLEGACY_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBLEGACY_SRC)))

LIBLEGACY_DEP := \
		$(LIBLEGACY_OBJ:.o=.d)

LIBLEGACY     := $(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBLEGACY)

clean-$(MODNAME)-dep:
		-rm -f $(LIBLEGACY_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBLEGACY_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(LIBLEGACY)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBLEGACY): $(LIBLEGACY_OBJ)
		$(MAKELIB) $@ $^

# add boost and eigen flags for the test object files and dependencies
$(LIBLEGACY_OBJ) $(LIBLEGACY_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ALLDEP += $(LIBLEGACY_DEP)
ALLLIB += $(LIBLEGACY)
