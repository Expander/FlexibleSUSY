DIR          := @DIR@
MODNAME      := lib@MODEL@

LIB@MODEL@_SRC := 

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIB@MODEL@_SRC += 
endif

ifneq ($(findstring lattice,$(ALGORITHMS)),)
LIB@MODEL@_SRC += 
endif

LIB@MODEL@_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB@MODEL@_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIB@MODEL@_SRC)))

LIB@MODEL@_DEP := \
		$(LIB@MODEL@_OBJ:.o=.d)

LIB@MODEL@     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIB@MODEL@)

clean-$(MODNAME):
		rm -rf $(LIB@MODEL@_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIB@MODEL@_DEP)
		rm -rf $(LIB@MODEL@)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LIB@MODEL@_DEP) $(LIB@MODEL@_OBJ): CPPFLAGS += $(TVMETFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)
endif

$(LIB@MODEL@): $(LIB@MODEL@_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIB@MODEL@_DEP)
ALLLIB += $(LIB@MODEL@)
