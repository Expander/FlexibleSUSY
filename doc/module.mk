DIR          := doc
MODNAME      := doc

HTML_OUTPUT_DIR := $(DIR)/html
PDF_OUTPUT_DIR  := $(DIR)
INDEX_PADE      := $(HTML_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile
DOXYGEN_MAINPAGE:= $(DIR)/mainpage.dox
MANUAL_PDF      := $(PDF_OUTPUT_DIR)/flexiblesusy.pdf
MANUAL_SRC      := \
		$(DIR)/flexiblesusy.tex \
		$(DIR)/version.tex \
		$(DIR)/chapters/overview.tex \
		$(DIR)/chapters/quick_start.tex \
		$(DIR)/chapters/usage.tex \
		$(DIR)/chapters/output.tex
PAPER_PDF       := $(PDF_OUTPUT_DIR)/paper.pdf
PAPER_SRC       := $(DIR)/paper.tex

LATEX_TMP       := \
		$(patsubst %.pdf, %.aux, $(MANUAL_PDF) $(PAPER_PDF)) \
		$(patsubst %.pdf, %.log, $(MANUAL_PDF) $(PAPER_PDF)) \
		$(patsubst %.pdf, %.toc, $(MANUAL_PDF) $(PAPER_PDF)) \
		$(patsubst %.pdf, %.out, $(MANUAL_PDF) $(PAPER_PDF)) \
		$(patsubst %.pdf, %.spl, $(MANUAL_PDF) $(PAPER_PDF))

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PADE) doc doc-html doc-pdf

doc: all-$(MODNAME)

doc-pdf: $(MANUAL_PDF) $(PAPER_PDF)

doc-html: $(INDEX_PADE)

all-$(MODNAME): doc-html doc-pdf

clean-$(MODNAME):
		-rm -f $(LATEX_TMP)

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -rf $(HTML_OUTPUT_DIR)
		-rm -f $(DOXYGEN_MAINPAGE)
		-rm -f $(MANUAL_PDF) $(PAPER_PDF)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PADE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES)" ; \
		  echo "OUTPUT_DIRECTORY = $(HTML_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)"; \
		  echo "EXCLUDE_PATTERNS = */test/*"; \
		) | doxygen -

$(MANUAL_PDF): $(MANUAL_SRC)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		bibtex $(shell echo $< | rev | cut -d. -f2 | rev)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<

$(PAPER_PDF): $(PAPER_SRC)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
