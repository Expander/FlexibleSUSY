DIR          := config
MODNAME      := config

$(DIR)/config.h: Makefile
	rm -f $@-t $@
	{ echo '/* DO NOT EDIT! GENERATED AUTOMATICALLY! */'; \
	  echo '#define VERSION "$(VERSION)/"'; \
	  echo '#define PKGNAME "$(PKGNAME)/"'; \
	} | sed '/""/d' > $@-t
	mv $@-t $@

CONFIG_HDR := \
		$(DIR)/config.h

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(CONFIG_HDR)

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(CONFIG_HDR)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ALLHDR += $(CONFIG_HDR)
