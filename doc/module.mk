DIR          := doc
MODNAME      := doc

DOC_OUTPUT_DIR  := $(DIR)/html
INDEX_PADE      := $(DOC_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile
MANUAL_PDF      := $(DIR)/flexiblesusy.pdf
MANUAL_SRC      := \
		$(DIR)/flexiblesusy.tex \
		$(DIR)/chapters/overview.tex \
		$(DIR)/chapters/quick_start.tex \
		$(DIR)/chapters/usage.tex \
		$(DIR)/chapters/output.tex

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PADE) doc

doc: all-$(MODNAME)

doc-pdf: $(MANUAL_PDF)

all-$(MODNAME): $(INDEX_PADE) $(MANUAL_PDF)

clean-$(MODNAME):
		rm -f $(DIR)/*.aux
		rm -f $(DIR)/*.log
		rm -f $(DIR)/*.toc

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(DOC_OUTPUT_DIR)
		rm -f $(MANUAL_PDF)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PADE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES)" ; \
		  echo "OUTPUT_DIRECTORY = $(DOC_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)" \
		) | doxygen -

$(MANUAL_PDF): $(MANUAL_SRC)
		pdflatex -output-directory $(DIR) $<
		pdflatex -output-directory $(DIR) $<
