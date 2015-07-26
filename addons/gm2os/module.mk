DIR          := addons/gm2os
MODNAME      := gm2os

LIBgm2os_INSTALL_DIR := gm2os

LIBgm2os_MK  := \
		$(DIR)/Makefile.in

LIBgm2os_gm2_MK  := \
		$(DIR)/module.gm2.mk

ifeq ($(shell $(FSCONFIG) --with-MSSMNoFVSLHA2),yes)
# source files
LIBgm2os_SRC := \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp

# main()
EXEgm2os_SRC := \
		$(DIR)/gm2.cpp

# header files
LIBgm2os_HDR := \
		$(DIR)/ffunctions.hpp \
		$(DIR)/gm2_1loop.hpp \
		$(DIR)/gm2_2loop.hpp \
		$(DIR)/gm2_error.hpp \
		$(DIR)/gm2_slha_io.hpp \
		$(DIR)/MSSMNoFV_onshell.hpp
endif

LIBgm2os_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBgm2os_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBgm2os_SRC)))

EXEgm2os_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEgm2os_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEgm2os_SRC)))

LIBgm2os_DEP := \
		$(LIBgm2os_OBJ:.o=.d)

EXEgm2os_DEP := \
		$(EXEgm2os_OBJ:.o=.d)

EXEgm2os_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEgm2os_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEgm2os_SRC)))

LIBgm2os     := \
		$(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         clean-$(MODNAME) clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) export-gm2os

clean-$(MODNAME)-dep:
		-rm -f $(LIBgm2os_DEP)
		-rm -f $(EXEgm2os_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBgm2os_OBJ)
		-rm -f $(EXEgm2os_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBgm2os)
		-rm -f $(EXEgm2os_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

export-gm2os:
		install -d $(LIBgm2os_INSTALL_DIR)
		install -d $(LIBgm2os_INSTALL_DIR)/gm2
		install -m u=rw,g=r,o=r $(LIBgm2os_SRC) $(LIBgm2os_INSTALL_DIR)/gm2
		install -m u=rw,g=r,o=r $(LIBgm2os_HDR) $(LIBgm2os_INSTALL_DIR)/gm2
		install -m u=rw,g=r,o=r $(LIBgm2os_MK) $(LIBgm2os_INSTALL_DIR)/Makefile
		install -m u=rw,g=r,o=r $(LIBgm2os_gm2_MK) $(LIBgm2os_INSTALL_DIR)/gm2/module.mk

$(LIBgm2os_DEP) $(EXEgm2os_DEP) $(LIBgm2os_OBJ) $(EXEgm2os_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBgm2os_DEP) $(EXEgm2os_DEP) $(LIBgm2os_OBJ) $(EXEgm2os_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBgm2os): $(LIBgm2os_OBJ)
		$(MAKELIB) $@ $^

$(EXEgm2os_EXE): $(EXEgm2os_OBJ) $(LIBgm2os) $(LIBMSSMNoFVSLHA2) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBgm2os_DEP) $(EXEgm2os_DEP)
ALLSRC += $(LIBgm2os_SRC) $(EXEgm2os_SRC)
ALLLIB += $(LIBgm2os)
ALLEXE += $(EXEgm2os_EXE)
