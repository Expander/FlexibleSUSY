DIR          := config
MODNAME      := config

CONFIG_HDR := \
		$(DIR)/config.h

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(CONFIG_HDR)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
