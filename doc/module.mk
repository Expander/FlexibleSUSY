DIR          := doc
MODNAME      := doc

DOC_MK       := \
		$(DIR)/module.mk

DOC_TMPL     := \
		$(DIR)/FlexibleEFTHiggs.rst \
		$(DIR)/HSSUSY.rst \
		$(DIR)/librarylink.rst \
		$(DIR)/meta_code.rst \
		$(DIR)/mainpage.dox.in \
		$(DIR)/model_file.rst \
		$(DIR)/MSSMEFTHiggs.rst \
		$(DIR)/NUHMSSMNoFVHimalaya.rst \
		$(DIR)/slha_input.rst

DOC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HTML_OUTPUT_DIR := $(DIR)/html
MAN_OUTPUT_DIR  := $(DIR)/man
PDF_OUTPUT_DIR  := $(DIR)
IMAGE_DIR       := $(DIR)/images
IMAGES          := $(IMAGE_DIR)/FS-logo.png \
		   $(IMAGE_DIR)/HSSUSY_Mh_MS.png \
		   $(IMAGE_DIR)/HSSUSY_Mh_Xt.png \
		   $(IMAGE_DIR)/HSSUSY_tower.svg \
		   $(IMAGE_DIR)/MSSMEFTHiggs_Mh_MS.png \
		   $(IMAGE_DIR)/MSSMEFTHiggs_tower.svg \
		   $(IMAGE_DIR)/NUHMSSMNoFVHimalaya_tower.svg
EXAMPLES_DIR    := $(DIR)/examples
EXAMPLES        := $(EXAMPLES_DIR)/HSSUSY_uncertainty_estimate.m \
		   $(EXAMPLES_DIR)/MSSMEFTHiggs_uncertainty_estimate.m \
		   $(EXAMPLES_DIR)/NUHMSSMNoFVHimalaya_uncertainty_estimate.m
INDEX_PAGE      := $(HTML_OUTPUT_DIR)/index.html
MAN_PAGE        := $(MAN_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile
DOXYGEN_MAINPAGE:= $(DIR)/mainpage.dox

PAPER_PDF_1     := $(PDF_OUTPUT_DIR)/flexiblesusy-1.0.pdf
PAPER_PDF_2     := $(PDF_OUTPUT_DIR)/flexiblesusy-2.0.pdf
PAPER_PDF       := $(PAPER_PDF_1) $(PAPER_PDF_2)
PAPER_SRC_1     := $(DIR)/flexiblesusy-1.0.tex
PAPER_SRC_2     := $(DIR)/flexiblesusy-2.0.tex
PAPER_SRC       := $(PAPER_SRC_1) $(PAPER_SRC_2)
PAPER_STY       := $(DIR)/JHEP.bst $(DIR)/tikz-uml.sty

LATEX_TMP       := \
		$(patsubst %.pdf, %.aux, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.bbl, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.blg, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.log, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.toc, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.out, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.spl, $(PAPER_PDF))

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PAGE) $(MAN_PAGE) doc doc-html doc-man doc-pdf

doc: all-$(MODNAME)

doc-pdf: $(PAPER_PDF)

doc-html: $(INDEX_PAGE)

doc-man: $(MAN_PAGE)

all-$(MODNAME): doc-html doc-man doc-pdf
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(DOC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(DOC_TMPL) $(DOC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(DOC_MK) $(DOC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(DOXYFILE) $(DOC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(PAPER_SRC) $(DOC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(PAPER_STY) $(DOC_INSTALL_DIR)
		install -d $(INSTALL_DIR)/$(IMAGE_DIR)
		install -m u=rw,g=r,o=r $(IMAGES) $(INSTALL_DIR)/$(IMAGE_DIR)
		install -d $(INSTALL_DIR)/$(EXAMPLES_DIR)
		install -m u=rw,g=r,o=r $(EXAMPLES) $(INSTALL_DIR)/$(EXAMPLES_DIR)
endif

clean-$(MODNAME):
		-rm -f $(LATEX_TMP)

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -rf $(HTML_OUTPUT_DIR)
		-rm -f $(DOXYGEN_MAINPAGE)
		-rm -f $(PAPER_PDF)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PAGE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES) $(README_FILE)" ; \
		  echo "OUTPUT_DIRECTORY = $(HTML_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)"; \
		  echo "EXCLUDE_PATTERNS = */meta/* */test/*"; \
		  echo "IMAGE_PATH = $(IMAGE_DIR)"; \
		  echo "INCLUDE_PATH = $(MODULES)"; \
		  echo "EXAMPLE_PATH = $(EXAMPLES_DIR)"; \
		) | doxygen -

$(MAN_PAGE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES) $(README_FILE)" ; \
		  echo "OUTPUT_DIRECTORY = $(MAN_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)"; \
		  echo "EXCLUDE_PATTERNS = */meta/* */test/*"; \
		  echo "IMAGE_PATH = $(IMAGE_DIR)"; \
		  echo "INCLUDE_PATH = $(MODULES)"; \
		  echo "EXAMPLE_PATH = $(EXAMPLES_DIR)"; \
		  echo "GENERATE_MAN = YES"; \
		  echo "GENERATE_HTML = NO"; \
		) | doxygen -

$(PAPER_PDF_1): $(PAPER_SRC_1) $(PAPER_STY)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<

$(PAPER_PDF_2): $(PAPER_SRC_2) $(PAPER_STY)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		bibtex $(PAPER_SRC_2:.tex=)
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
		pdflatex -output-directory $(PDF_OUTPUT_DIR) $<
