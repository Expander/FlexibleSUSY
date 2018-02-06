DIR          := higher_order/SM
MODNAME      := higher_order_SM
WITH_$(MODNAME) := yes

LIB_higher_order_SM_MK  := \
		$(DIR)/module.mk

LIB_higher_order_SM_SRC := \
		$(DIR)/sm_fourloophiggs.cpp \
		$(DIR)/sm_threeloop_as.cpp \
		$(DIR)/sm_threeloophiggs.cpp \
		$(DIR)/sm_twoloophiggs.cpp

LIB_higher_order_SM_HDR := \
		$(DIR)/sm_fourloophiggs.hpp \
		$(DIR)/sm_threeloop_as.hpp \
		$(DIR)/sm_threeloophiggs.hpp \
		$(DIR)/sm_twoloophiggs.hpp

LIB_higher_order_SM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_higher_order_SM_SRC)))

LIB_higher_order_SM_DEP := $(LIB_higher_order_SM_OBJ:.o=.d)

LIB_higher_order_SM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_higher_order_SM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_higher_order_SM)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIB_higher_order_SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SM_SRC) $(LIB_higher_order_SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SM_HDR) $(LIB_higher_order_SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_higher_order_SM_MK) $(LIB_higher_order_SM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIB_higher_order_SM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIB_higher_order_SM)

clean-$(MODNAME)-obj:
		-rm -f $(LIB_higher_order_SM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_higher_order_SM_DEP) $(LIB_higher_order_SM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_higher_order_SM_DEP) $(LIB_higher_order_SM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_higher_order_SM): $(LIB_higher_order_SM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIB_higher_order_SM): $(LIB_higher_order_SM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_higher_order_SM_DEP)
ALLLIB += $(LIB_higher_order_SM)
