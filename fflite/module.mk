DIR          := fflite
MODNAME      := libfflite

LIBFFLITE_SRC := \
		$(DIR)/ffca0.F \
		$(DIR)/ffcb0.F \
		$(DIR)/ffcb1.F \
		$(DIR)/ffcb2p.F \
		$(DIR)/ffcc0.F \
		$(DIR)/ffcli2.F \
		$(DIR)/ffinit.F \
		$(DIR)/ffxa0.F \
		$(DIR)/ffxb0.F \
		$(DIR)/ffxb1.F \
		$(DIR)/ffxb2p.F \
		$(DIR)/ffxli2.F \
		$(DIR)/ini.F

LIBFFLITE_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFFLITE_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFFLITE_SRC))) \
		$(patsubst %.F, %.o, $(filter %.F, $(LIBFFLITE_SRC)))

LIBFFLITE_DEP := \
		$(LIBFFLITE_OBJ:.o=.d)

LIBFFLITE     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFFLITE)

clean-$(MODNAME):
		rm -rf $(LIBFFLITE_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFFLITE_DEP)
		rm -rf $(LIBFFLITE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFFLITE): $(LIBFFLITE_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFFLITE_DEP)
ALLLIB += $(LIBFFLITE)
