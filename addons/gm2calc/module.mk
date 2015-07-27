DIR          := addons/gm2calc
MODNAME      := gm2calc

LIBgm2calc_INSTALL_DIR := gm2calc

LIBgm2calc_MK  := \
		$(DIR)/Makefile.in

LIBgm2calc_gm2_MK  := \
		$(DIR)/module.src.mk

# source files
LIBgm2calc_SRC := \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFV_onshell_physical.cpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.cpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.cpp

# main()
EXEgm2calc_SRC := \
		$(DIR)/gm2.cpp

# header files
LIBgm2calc_HDR := \
		$(DIR)/ffunctions.hpp \
		$(DIR)/gm2_1loop.hpp \
		$(DIR)/gm2_2loop.hpp \
		$(DIR)/gm2_error.hpp \
		$(DIR)/gm2_slha_io.hpp \
		$(DIR)/MSSMNoFV_onshell.hpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFV_onshell_physical.hpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.hpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.hpp

LIBgm2calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBgm2calc_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBgm2calc_SRC)))

EXEgm2calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEgm2calc_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEgm2calc_SRC)))

LIBgm2calc_DEP := \
		$(LIBgm2calc_OBJ:.o=.d)

EXEgm2calc_DEP := \
		$(EXEgm2calc_OBJ:.o=.d)

EXEgm2calc_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEgm2calc_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEgm2calc_SRC)))

LIBgm2calc     := \
		$(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         clean-$(MODNAME) clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) export-gm2calc

clean-$(MODNAME)-dep:
		-rm -f $(LIBgm2calc_DEP)
		-rm -f $(EXEgm2calc_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBgm2calc_OBJ)
		-rm -f $(EXEgm2calc_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBgm2calc)
		-rm -f $(EXEgm2calc_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

# dependent files in src/
LIBgm2calc_DEP_SRC := \
		src/compare.hpp \
		src/dilog.h \
		src/dilog.f \
		src/eigen_utils.hpp \
		src/error.hpp \
		src/linalg2.hpp \
		src/numerics2.hpp

export-gm2calc:
		install -d $(LIBgm2calc_INSTALL_DIR)
		install -d $(LIBgm2calc_INSTALL_DIR)/src
		install -m u=rw,g=r,o=r $(LIBgm2calc_MK) $(LIBgm2calc_INSTALL_DIR)/Makefile
		install -m u=rw,g=r,o=r $(LIBgm2calc_gm2_MK) $(LIBgm2calc_INSTALL_DIR)/src/module.mk
		install -m u=rw,g=r,o=r $(LIBgm2calc_SRC) $(LIBgm2calc_INSTALL_DIR)/src
		install -m u=rw,g=r,o=r $(LIBgm2calc_HDR) $(LIBgm2calc_INSTALL_DIR)/src
		install -m u=rw,g=r,o=r $(EXEgm2calc_SRC) $(LIBgm2calc_INSTALL_DIR)/src
		install -m u=rw,g=r,o=r $(LIBgm2calc_DEP_SRC) $(LIBgm2calc_INSTALL_DIR)/src
		install -m u=rw,g=r,o=r slhaea/slhaea.h $(LIBgm2calc_INSTALL_DIR)/src

$(LIBgm2calc_DEP) $(EXEgm2calc_DEP) $(LIBgm2calc_OBJ) $(EXEgm2calc_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(LIBgm2calc): $(LIBgm2calc_OBJ)
		$(MAKELIB) $@ $^

$(EXEgm2calc_EXE): $(EXEgm2calc_OBJ) $(LIBgm2calc) $(LIBMSSMNoFV_onshell) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBgm2calc_DEP) $(EXEgm2calc_DEP)
ALLSRC += $(LIBgm2calc_SRC) $(EXEgm2calc_SRC)
ALLLIB += $(LIBgm2calc)
ALLEXE += $(EXEgm2calc_EXE)
