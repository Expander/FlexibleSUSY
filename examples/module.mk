DIR      := examples
MODNAME  := examples

EXAMPLES_SRC := \
		$(DIR)/softsusy.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
EXAMPLES_SRC += \
		$(DIR)/run_mssm.cpp
endif

EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXAMPLES_SRC)))

EXAMPLES_DEP := \
		$(EXAMPLES_OBJ:.o=.d)

EXAMPLES_EXE := \
		$(EXAMPLES_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(EXAMPLES_EXE)

clean-$(MODNAME):
		rm -rf $(EXAMPLES_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(EXAMPLES_DEP)
		rm -rf $(EXAMPLES_EXE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/run_mssm.x: $(DIR)/run_mssm.o $(LIBMSSM) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(FLIBS)

$(DIR)/softsusy.x: $(DIR)/softsusy.o $(LIBFLEXI) $(LIBMSSM)
		$(CXX) -o $@ $^ $(FLIBS)

ALLDEP += $(EXAMPLES_DEP)
ALLEXE += $(EXAMPLES_EXE)
