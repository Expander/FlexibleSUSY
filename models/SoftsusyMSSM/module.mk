DIR          := models/SoftsusyMSSM
MODNAME      := SoftsusyMSSM
WITH_$(MODNAME) := yes
MODSoftsusyMSSM_MOD := NMSSM_higgs MSSM_higgs
MODSoftsusyMSSM_DEP := $(patsubst %,model_specific/%,$(MODSoftsusyMSSM_MOD))
MODSoftsusyMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSoftsusyMSSM_MOD))
MODSoftsusyMSSM_LIB := $(foreach M,$(MODSoftsusyMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

LIBSoftsusyMSSM_HDR  := \
		$(DIR)/conversion.hpp \
		$(DIR)/def.h \
		$(DIR)/diagonalization.hpp \
		$(DIR)/gut_scale_calculator.hpp \
		$(DIR)/linalg.h \
		$(DIR)/lowe_legacy.h \
		$(DIR)/mycomplex.h \
		$(DIR)/numerics_legacy.h \
		$(DIR)/rge.h \
		$(DIR)/rk_legacy.hpp \
		$(DIR)/utils.h \
		$(DIR)/xpr-base.h \
		$(DIR)/xpr-matrix.h \
		$(DIR)/xpr-vector.h

LIBSoftsusyMSSM_SRC  := \
		$(DIR)/conversion.cpp \
		$(DIR)/def.cpp \
		$(DIR)/diagonalization.cpp \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe_legacy.cpp \
		$(DIR)/mssmUtils.cpp \
		$(DIR)/numerics_legacy.cpp \
		$(DIR)/physpars.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/rk_legacy.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/tensor.cpp \
		$(DIR)/utils.cpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
LIBSoftsusyMSSM_SRC  += \
		$(DIR)/SoftsusyMSSM_two_scale.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_susy_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_sugra_constraint.cpp
endif

LIBSoftsusyMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSoftsusyMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSoftsusyMSSM_SRC)))

LIBSoftsusyMSSM_DEP  := \
		$(LIBSoftsusyMSSM_OBJ:.o=.d)

LIBSoftsusyMSSM      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

EXESoftsusyMSSM_SRC := \
		$(DIR)/run_SoftsusyMSSM.cpp

EXESoftsusyMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESoftsusyMSSM_SRC)))

EXESoftsusyMSSM_DEP := \
		$(EXESoftsusyMSSM_OBJ:.o=.d)

RUN_SoftsusyMSSM_EXE := \
		$(EXESoftsusyMSSM_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj

all-$(MODNAME): $(LIBSoftsusyMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBSoftsusyMSSM_DEP)
		-rm -f $(EXESoftsusyMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSoftsusyMSSM_OBJ)
		-rm -f $(EXESoftsusyMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSoftsusyMSSM)
		-rm -f $(RUN_SoftsusyMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSoftsusyMSSM_DEP) $(EXESoftsusyMSSM_DEP) $(LIBSoftsusyMSSM_OBJ) $(EXESoftsusyMSSM_OBJ): \
	CPPFLAGS += $(MODSoftsusyMSSM_INC) $(EIGENFLAGS)

$(LIBSoftsusyMSSM): $(LIBSoftsusyMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(RUN_SoftsusyMSSM_EXE): $(EXESoftsusyMSSM_OBJ) $(LIBSoftsusyMSSM) $(MODSoftsusyMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

ALLDEP += $(LIBSoftsusyMSSM_DEP) $(EXESoftsusyMSSM_DEP)
ALLLIB += $(LIBSoftsusyMSSM)
ALLEXE += $(RUN_SoftsusyMSSM_EXE)
ALLMODDEP += $(MODSoftsusyMSSM_DEP)
