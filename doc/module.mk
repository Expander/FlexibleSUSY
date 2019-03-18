DIR          := doc
MODNAME      := doc

DOC_MK       := \
		$(DIR)/module.mk

MODELS_DIR    :=$(DIR)/models
DOC_TMPL     := \
		$(MODELS_DIR)/FlexibleEFTHiggs.rst \
		$(MODELS_DIR)/HSSUSY.rst \
		$(MODELS_DIR)/librarylink.rst \
		$(MODELS_DIR)/meta_code.rst \
		$(MODELS_DIR)/mainpage.dox.in \
		$(MODELS_DIR)/model_file.rst \
		$(MODELS_DIR)/MSSMEFTHiggs.rst \
		$(MODELS_DIR)/NUHMSSMNoFVHimalaya.rst \
		$(MODELS_DIR)/slha_input.rst

DOC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HTML_OUTPUT_DIR := $(DIR)/html
MAN_OUTPUT_DIR  := $(DIR)/man
MANUALS_DIR     := $(DIR)/manuals
PDF_OUTPUT_DIR  := $(MANUALS_DIR)


IMAGE_DIR       := $(MANUALS_DIR)/images
IMAGES          := $(IMAGE_DIR)/FS-logo.png \
		   $(IMAGE_DIR)/HSSUSY_Mh_MS.png \
		   $(IMAGE_DIR)/HSSUSY_Mh_Xt.png \
		   $(IMAGE_DIR)/HSSUSY_tower.svg \
		   $(IMAGE_DIR)/MSSMEFTHiggs_Mh_MS.png \
		   $(IMAGE_DIR)/MSSMEFTHiggs_tower.svg \
		   $(IMAGE_DIR)/NUHMSSMNoFVHimalaya_tower.svg
EXAMPLES_DIR    := $(MANUALS_DIR)/examples
EXAMPLES        := $(EXAMPLES_DIR)/HSSUSY_uncertainty_estimate.m \
		   $(EXAMPLES_DIR)/MSSMEFTHiggs_uncertainty_estimate.m \
		   $(EXAMPLES_DIR)/NUHMSSMNoFVHimalaya_uncertainty_estimate.m
INDEX_PAGE      := $(HTML_OUTPUT_DIR)/index.html
MAN_PAGE        := $(MAN_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile
DOXYGEN_MAINPAGE:= $(DIR)/mainpage.dox

PAPER_PDF_1     := $(PDF_OUTPUT_DIR)/flexiblesusy-1.0.pdf
PAPER_PDF_2     := $(PDF_OUTPUT_DIR)/flexiblesusy-2.0.pdf
PAPER_PDF_3     := $(PDF_OUTPUT_DIR)/flexiblesusy-new_features.pdf
PAPER_PDF       := $(PAPER_PDF_1) $(PAPER_PDF_2) $(PAPER_PDF_3)
PAPER_SRC_1     := $(MANUALS_DIR)/flexiblesusy-1.0.tex
PAPER_SRC_2     := $(MANUALS_DIR)/flexiblesusy-2.0.tex
PAPER_SRC_3     := $(MANUALS_DIR)/flexiblesusy-new_features.tex
PAPER_SRC       := $(PAPER_SRC_1) $(PAPER_SRC_2) $(PAPER_SRC_3)
PAPER_STY       := $(MANUALS_DIR)/JHEP.bst $(MANUALS_DIR)/tikz-uml.sty

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
		cd doc/manuals && pdflatex flexiblesusy-1.0.tex
		cd doc/manuals && pdflatex flexiblesusy-1.0.tex
		cd doc/manuals && pdflatex flexiblesusy-1.0.tex

$(PAPER_PDF_2): $(PAPER_SRC_2) $(PAPER_STY)
		cd doc/manuals && pdflatex flexiblesusy-2.0.tex
		cd doc/manuals && bibtex flexiblesusy-2.0
		cd doc/manuals && pdflatex flexiblesusy-2.0.tex
		cd doc/manuals && pdflatex flexiblesusy-2.0.tex

$(PAPER_PDF_3): $(PAPER_SRC_3) $(PAPER_STY)
		cd doc/manuals && pdflatex flexiblesusy-new_features.tex
		cd doc/manuals && bibtex flexiblesusy-new_features
		cd doc/manuals && pdflatex flexiblesusy-new_features.tex
		cd doc/manuals && pdflatex flexiblesusy-new_features.tex
