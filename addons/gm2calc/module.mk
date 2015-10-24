DIR          := addons/gm2calc
MODNAME      := gm2calc

# source files
LIBgm2calc_SRC := \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_mb.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/gm2_uncertainty.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFV_onshell_physical.cpp \
		$(DIR)/MSSMNoFV_onshell_problems.cpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.cpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.cpp

# main()
EXEgm2calc_SRC := \
		$(DIR)/gm2calc.cpp \
		$(DIR)/gm2scan.cpp

# header files
LIBgm2calc_HDR := \
		$(DIR)/ffunctions.hpp \
		$(DIR)/gm2_1loop.hpp \
		$(DIR)/gm2_2loop.hpp \
		$(DIR)/gm2_error.hpp \
		$(DIR)/gm2_mb.hpp \
		$(DIR)/gm2_slha_io.hpp \
		$(DIR)/gm2_uncertainty.hpp \
		$(DIR)/MSSMNoFV_onshell.hpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFV_onshell_physical.hpp \
		$(DIR)/MSSMNoFV_onshell_problems.hpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.hpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.hpp

LIBgm2calc_SLHA_INPUT := \
		$(DIR)/example.slha \
		$(DIR)/example.gm2

LIBgm2calc_INFO := \
		$(DIR)/AUTHORS \
		$(DIR)/COPYING \
		$(DIR)/README

LIBgm2calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBgm2calc_SRC)))

EXEgm2calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEgm2calc_SRC)))

LIBgm2calc_DEP := \
		$(LIBgm2calc_OBJ:.o=.d)

EXEgm2calc_DEP := \
		$(EXEgm2calc_OBJ:.o=.d)

EXEgm2calc_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEgm2calc_SRC)))

LIBgm2calc     := \
		$(DIR)/lib$(MODNAME)$(LIBEXT)

.PHONY:         clean-$(MODNAME) clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME)

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

$(LIBgm2calc_DEP) $(EXEgm2calc_DEP) $(LIBgm2calc_OBJ) $(EXEgm2calc_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(LIBgm2calc): $(LIBgm2calc_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBgm2calc) $(LIBFLEXI)
		$(CXX) -o $@ $(call abspathx,$^)

ALLDEP += $(LIBgm2calc_DEP) $(EXEgm2calc_DEP)
ALLSRC += $(LIBgm2calc_SRC) $(EXEgm2calc_SRC)
ALLLIB += $(LIBgm2calc)
ALLEXE += $(EXEgm2calc_EXE)
