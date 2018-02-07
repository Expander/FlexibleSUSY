DIR          := higher_order/MSSM_higgs
MODNAME      := higher_order_MSSM_higgs
WITH_$(MODNAME) := yes

LIB_higher_order_MSSM_higgs_MK  := \
		$(DIR)/module.mk

LIB_higher_order_MSSM_higgs_SRC := \
		$(DIR)/mssm_twoloophiggs.cpp \
		$(DIR)/mssm_twoloophiggs_impl.f

LIB_higher_order_MSSM_higgs_HDR := \
		$(DIR)/mssm_twoloophiggs.h \
		$(DIR)/mssm_twoloophiggs.hpp

LIB_higher_order_MSSM_higgs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_higher_order_MSSM_higgs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIB_higher_order_MSSM_higgs_SRC)))

LIB_higher_order_MSSM_higgs_DEP := $(LIB_higher_order_MSSM_higgs_OBJ:.o=.d)

LIB_higher_order_MSSM_higgs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_higher_order_MSSM_higgs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_higher_order_MSSM_higgs)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIB_higher_order_MSSM_higgs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_higgs_SRC) $(LIB_higher_order_MSSM_higgs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_higgs_HDR) $(LIB_higher_order_MSSM_higgs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_higgs_MK) $(LIB_higher_order_MSSM_higgs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIB_higher_order_MSSM_higgs_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIB_higher_order_MSSM_higgs)

clean-$(MODNAME)-obj:
		-rm -f $(LIB_higher_order_MSSM_higgs_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_higher_order_MSSM_higgs_DEP) $(LIB_higher_order_MSSM_higgs_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_higher_order_MSSM_higgs_DEP) $(LIB_higher_order_MSSM_higgs_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_higher_order_MSSM_higgs): $(LIB_higher_order_MSSM_higgs_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIB_higher_order_MSSM_higgs): $(LIB_higher_order_MSSM_higgs_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_higher_order_MSSM_higgs_DEP)
ALLLIB += $(LIB_higher_order_MSSM_higgs)
