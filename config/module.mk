DIR          := config
MODNAME      := config

CONFIG_HDR := \
		$(DIR)/config.h

REQUIRED_SARAH_VERSION_FILE := \
		$(DIR)/required_sarah_version.m

FLEXIBLESUSY_VERSION_FILE := \
		$(DIR)/version

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(CONFIG_HDR)
		-rm -f $(FLEXIBLESUSY_VERSION_FILE)
		-rm -f $(REQUIRED_SARAH_VERSION_FILE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
