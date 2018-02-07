DIR          := higher_order/MSSM_thresholds
MODNAME      := higher_order_MSSM_thresholds
WITH_$(MODNAME) := yes

LIB_higher_order_MSSM_thresholds_MK  := \
		$(DIR)/module.mk

LIB_higher_order_MSSM_thresholds_SRC := \
		$(DIR)/mssm_twoloop_as.cpp \
		$(DIR)/mssm_twoloop_mb.cpp \
		$(DIR)/mssm_twoloop_mt.cpp \
		$(DIR)/mssm_twoloop_mtau.cpp

LIB_higher_order_MSSM_thresholds_HDR := \
		$(DIR)/mssm_twoloop_as.hpp \
		$(DIR)/mssm_twoloop_mb.hpp \
		$(DIR)/mssm_twoloop_mt.hpp \
		$(DIR)/mssm_twoloop_mtau.hpp

LIB_higher_order_MSSM_thresholds_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_higher_order_MSSM_thresholds_SRC)))

LIB_higher_order_MSSM_thresholds_DEP := $(LIB_higher_order_MSSM_thresholds_OBJ:.o=.d)

LIB_higher_order_MSSM_thresholds     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_higher_order_MSSM_thresholds_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_higher_order_MSSM_thresholds)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIB_higher_order_MSSM_thresholds_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_thresholds_SRC) $(LIB_higher_order_MSSM_thresholds_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_thresholds_HDR) $(LIB_higher_order_MSSM_thresholds_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_MSSM_thresholds_MK) $(LIB_higher_order_MSSM_thresholds_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIB_higher_order_MSSM_thresholds_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIB_higher_order_MSSM_thresholds)

clean-$(MODNAME)-obj:
		-rm -f $(LIB_higher_order_MSSM_thresholds_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_higher_order_MSSM_thresholds_DEP) $(LIB_higher_order_MSSM_thresholds_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_higher_order_MSSM_thresholds_DEP) $(LIB_higher_order_MSSM_thresholds_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_higher_order_MSSM_thresholds): $(LIB_higher_order_MSSM_thresholds_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIB_higher_order_MSSM_thresholds): $(LIB_higher_order_MSSM_thresholds_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_higher_order_MSSM_thresholds_DEP)
ALLLIB += $(LIB_higher_order_MSSM_thresholds)
