DIR          := H3m
MODNAME      := H3m
WITH_$(MODNAME) := yes

LIBH3M_INC := \
		$(DIR)/hierarchies/h32q2g/sigS12Full.inc \
		$(DIR)/hierarchies/h32q2g/sigS1Full.inc \
		$(DIR)/hierarchies/h32q2g/sigS2Full.inc \
		$(DIR)/hierarchies/h3q22g/sigS12Full.inc \
		$(DIR)/hierarchies/h3q22g/sigS1Full.inc \
		$(DIR)/hierarchies/h3q22g/sigS2Full.inc \
		$(DIR)/hierarchies/h3/sigS12Full.inc \
		$(DIR)/hierarchies/h3/sigS1Full.inc \
		$(DIR)/hierarchies/h3/sigS2Full.inc \
		$(DIR)/hierarchies/h4/sigS12Full.inc \
		$(DIR)/hierarchies/h4/sigS1Full.inc \
		$(DIR)/hierarchies/h4/sigS2Full.inc \
		$(DIR)/hierarchies/h5g1/sigS12Full.inc \
		$(DIR)/hierarchies/h5g1/sigS1Full.inc \
		$(DIR)/hierarchies/h5g1/sigS2Full.inc \
		$(DIR)/hierarchies/h5/sigS12Full.inc \
		$(DIR)/hierarchies/h5/sigS1Full.inc \
		$(DIR)/hierarchies/h5/sigS2Full.inc \
		$(DIR)/hierarchies/h6b2qg2/sigS12Full.inc \
		$(DIR)/hierarchies/h6b2qg2/sigS1Full.inc \
		$(DIR)/hierarchies/h6b2qg2/sigS2Full.inc \
		$(DIR)/hierarchies/h6bq22g/sigS12Full.inc \
		$(DIR)/hierarchies/h6bq22g/sigS1Full.inc \
		$(DIR)/hierarchies/h6bq22g/sigS2Full.inc \
		$(DIR)/hierarchies/h6bq2g2/sigS12Full.inc \
		$(DIR)/hierarchies/h6bq2g2/sigS1Full.inc \
		$(DIR)/hierarchies/h6bq2g2/sigS2Full.inc \
		$(DIR)/hierarchies/h6b/sigS12Full.inc \
		$(DIR)/hierarchies/h6b/sigS1Full.inc \
		$(DIR)/hierarchies/h6b/sigS2Full.inc \
		$(DIR)/hierarchies/h6g2/sigS12Full.inc \
		$(DIR)/hierarchies/h6g2/sigS1Full.inc \
		$(DIR)/hierarchies/h6g2/sigS2Full.inc \
		$(DIR)/hierarchies/h6/sigS12Full.inc \
		$(DIR)/hierarchies/h6/sigS1Full.inc \
		$(DIR)/hierarchies/h6/sigS2Full.inc \
		$(DIR)/hierarchies/h9q2/sigS12Full.inc \
		$(DIR)/hierarchies/h9q2/sigS1Full.inc \
		$(DIR)/hierarchies/h9q2/sigS2Full.inc \
		$(DIR)/hierarchies/h9/sigS12Full.inc \
		$(DIR)/hierarchies/h9/sigS1Full.inc \
		$(DIR)/hierarchies/h9/sigS2Full.inc \

LIBH3M_HDR := \
		$(DIR)/H3m_interface.hpp \
		$(DIR)/HierarchyCalculator.hpp \
		$(LIBH3M_INC)

LIBH3M_MK := \
		$(DIR)/module.mk

LIBH3M_SRC := \
		$(DIR)/HierarchyCalculator.cpp

LIBH3M_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBH3M_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBH3M_SRC))) \
		$(patsubst %.F, %.o, $(filter %.F, $(LIBH3M_SRC)))

LIBH3M_DEP := \
		$(LIBH3M_OBJ:.o=.d)

LIBH3M     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIBH3M_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBH3M)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBH3M_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBH3M_SRC) $(LIBH3M_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBH3M_HDR) $(LIBH3M_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBH3M_MK) $(LIBH3M_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBH3M_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBH3M)

clean-$(MODNAME)-obj:
		-rm -f $(LIBH3M_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBH3M_DEP) $(LIBH3M_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBH3M): $(LIBH3M_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBH3M_DEP)
ALLLIB += $(LIBH3M)
