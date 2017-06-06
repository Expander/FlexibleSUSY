DIR          := Himalaya
MODNAME      := Himalaya
WITH_$(MODNAME) := yes

LIBHIMALAYA_HIERARCHIES_DIR := $(DIR)/hierarchies

LIBHIMALAYA_INC := \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS2Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS12Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS1Full.inc \
		$(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS2Full.inc \

LIBHIMALAYA_HDR := \
		$(DIR)/Himalaya_interface.hpp \
		$(DIR)/HierarchyCalculator.hpp \
		$(DIR)/HierarchyObject.hpp

LIBHIMALAYA_MK := \
		$(DIR)/module.mk

LIBHIMALAYA_SRC := \
		$(DIR)/HierarchyCalculator.cpp \
		$(DIR)/HierarchyObject.cpp

LIBHIMALAYA_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBHIMALAYA_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBHIMALAYA_SRC))) \
		$(patsubst %.F, %.o, $(filter %.F, $(LIBHIMALAYA_SRC)))

LIBHIMALAYA_DEP := \
		$(LIBHIMALAYA_OBJ:.o=.d)

LIBHIMALAYA     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIBHIMALAYA_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBHIMALAYA)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBHIMALAYA_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_SRC) $(LIBHIMALAYA_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HDR) $(LIBHIMALAYA_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_MK) $(LIBHIMALAYA_INSTALL_DIR)
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h32q2g
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3q22g
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h4
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5g1
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b2qg2
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq22g
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq2g2
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6g2
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9
		install -d $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9q2
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS12Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h32q2g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS1Full.inc   $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h32q2g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h32q2g/sigS2Full.inc   $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h32q2g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS12Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3q22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS1Full.inc   $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3q22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3q22g/sigS2Full.inc   $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3q22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS12Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS1Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h3/sigS2Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h3/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS12Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h4/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS1Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h4/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h4/sigS2Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h4/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS12Full.inc    $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5g1/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS1Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5g1/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5g1/sigS2Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5g1/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS12Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS1Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h5/sigS2Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h5/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS12Full.inc $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b2qg2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS1Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b2qg2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b2qg2/sigS2Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b2qg2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS12Full.inc $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS1Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq22g/sigS2Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq22g/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS12Full.inc $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq2g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS1Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq2g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6bq2g2/sigS2Full.inc  $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6bq2g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS12Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS1Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6b/sigS2Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6b/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS12Full.inc    $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS1Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6g2/sigS2Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6g2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS12Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS1Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h6/sigS2Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h6/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS12Full.inc    $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9q2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS1Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9q2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9q2/sigS2Full.inc     $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9q2/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS12Full.inc      $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS1Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9/
		install -m u=rw,g=r,o=r $(LIBHIMALAYA_HIERARCHIES_DIR)/h9/sigS2Full.inc       $(LIBHIMALAYA_INSTALL_DIR)/hierarchies/h9/
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBHIMALAYA_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBHIMALAYA)

clean-$(MODNAME)-obj:
		-rm -f $(LIBHIMALAYA_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBHIMALAYA_DEP) $(LIBHIMALAYA_OBJ): CPPFLAGS += $(EIGENFLAGS)

$(LIBHIMALAYA): $(LIBHIMALAYA_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBHIMALAYA_DEP)
ALLLIB += $(LIBHIMALAYA)
