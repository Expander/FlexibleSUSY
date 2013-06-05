DIR          := doc
MODNAME      := doc

DOC_OUTPUT_DIR  := $(DIR)/html
INDEX_PADE      := $(DOC_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PADE) doc

doc: all-$(MODNAME)

all-$(MODNAME): $(INDEX_PADE)

clean-$(MODNAME):
		rm -rf $(DOC_OUTPUT_DIR)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PADE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES)" ; \
		  echo "OUTPUT_DIRECTORY = $(DOC_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP)" ) | doxygen -
