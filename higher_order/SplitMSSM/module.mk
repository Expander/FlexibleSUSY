DIR          := higher_order/SplitMSSM
MODNAME      := higher_order_SplitMSSM
WITH_$(MODNAME) := yes

LIB_higher_order_SplitMSSM_MK  := \
		$(DIR)/module.mk

LIB_higher_order_SplitMSSM_SRC := \
		$(DIR)/splitmssm_threeloophiggs.cpp \
		$(DIR)/splitmssm_thresholds.cpp

LIB_higher_order_SplitMSSM_HDR := \
		$(DIR)/splitmssm_threeloophiggs.hpp \
		$(DIR)/splitmssm_thresholds.hpp

LIB_higher_order_SplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_higher_order_SplitMSSM_SRC)))

LIB_higher_order_SplitMSSM_DEP := $(LIB_higher_order_SplitMSSM_OBJ:.o=.d)

LIB_higher_order_SplitMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_higher_order_SplitMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_higher_order_SplitMSSM)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIB_higher_order_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SplitMSSM_SRC) $(LIB_higher_order_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SplitMSSM_HDR) $(LIB_higher_order_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SplitMSSM_MK) $(LIB_higher_order_SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIB_higher_order_SplitMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIB_higher_order_SplitMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIB_higher_order_SplitMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_higher_order_SplitMSSM_DEP) $(LIB_higher_order_SplitMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_higher_order_SplitMSSM_DEP) $(LIB_higher_order_SplitMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_higher_order_SplitMSSM): $(LIB_higher_order_SplitMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIB_higher_order_SplitMSSM): $(LIB_higher_order_SplitMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_higher_order_SplitMSSM_DEP)
ALLLIB += $(LIB_higher_order_SplitMSSM)
