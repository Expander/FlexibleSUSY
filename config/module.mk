DIR          := config
MODNAME      := config

CONFIG_HDR := \
		$(DIR)/config.h

CONFIG_MK    := \
		$(DIR)/module.mk

CONFIG_TMPL  := \
		$(DIR)/abspathx.mk \
		$(DIR)/config.h.in \
		$(DIR)/flexiblesusy-config.in \
		$(DIR)/list_sarah_model_files.sh.in \
		$(DIR)/Makefile.in

CONFIG_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

REQUIRED_SARAH_VERSION_FILE := \
		$(DIR)/required_sarah_version.m

FLEXIBLESUSY_VERSION_FILE := \
		$(DIR)/version

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CONFIG_INSTALL_DIR)
		install $(CONFIG_TMPL) $(CONFIG_INSTALL_DIR)
		install $(CONFIG_MK) $(CONFIG_INSTALL_DIR)
endif

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(CONFIG_HDR)
		-rm -f $(FLEXIBLESUSY_VERSION_FILE)
		-rm -f $(REQUIRED_SARAH_VERSION_FILE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
