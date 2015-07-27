DIR          := gm2
MODNAME      := gm2

# source files
LIBgm2_SRC := \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp

# main()
EXEgm2_SRC := \
		$(DIR)/gm2.cpp

LIBgm2_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBgm2_SRC)))

EXEgm2_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEgm2_SRC)))

LIBgm2_DEP := \
		$(LIBgm2_OBJ:.o=.d)

EXEgm2_DEP := \
		$(EXEgm2_OBJ:.o=.d)

EXEgm2_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEgm2_SRC)))

LIBgm2     := \
		$(DIR)/lib$(MODNAME)$(LIBEXT)

clean::
		-rm -f $(LIBgm2_DEP)
		-rm -f $(EXEgm2_DEP)
		-rm -f $(LIBgm2_OBJ)
		-rm -f $(EXEgm2_OBJ)
		-rm -f $(LIBgm2)
		-rm -f $(EXEgm2_EXE)

$(LIBgm2_DEP) $(EXEgm2_DEP) $(LIBgm2_OBJ) $(EXEgm2_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(LIBgm2): $(LIBgm2_OBJ)
		$(MAKELIB) $@ $^

$(EXEgm2_EXE): $(EXEgm2_OBJ) $(LIBgm2) $(LIBMSSMNoFVSLHA2) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $^ $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBgm2_DEP) $(EXEgm2_DEP)
ALLLIB += $(LIBgm2)
ALLEXE += $(EXEgm2_EXE)
