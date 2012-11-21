DIR          := doc
MODNAME      := doc

INDEX_PADE      := $(DIR)/html/index.html
DOXYFILE        := $(DIR)/Doxyfile

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PADE)

all-$(MODNAME): $(INDEX_PADE)

clean-$(MODNAME):
		rm -rf $(DIR)/html/

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PADE):
		( cat $(DOXYFILE) ; \
                  echo "INPUT = $(MODULES)" ; \
                  echo "OUTPUT_DIRECTORY = $(DIR)/html" ) | doxygen -
