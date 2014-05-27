DIR      := examples
MODNAME  := examples

EXAMPLES_SRC :=

EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXAMPLES_SRC)))

EXAMPLES_DEP := \
		$(EXAMPLES_OBJ:.o=.d)

EXAMPLES_EXE := \
		$(EXAMPLES_OBJ:.o=.x)

STANDALONE_DIR := \
		$(DIR)/standalone-model \
		$(DIR)/standalone-rge

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(EXAMPLES_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(EXAMPLES_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(EXAMPLES_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(EXAMPLES_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		-@for d in $(STANDALONE_DIR); do \
			(cd $$d && make distclean); \
		 done

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ALLDEP += $(EXAMPLES_DEP)
ALLEXE += $(EXAMPLES_EXE)
