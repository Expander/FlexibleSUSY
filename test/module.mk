DIR      := test
MODNAME  := test

LIBTEST_SRC := \
		$(DIR)/stopwatch.cpp

LIBTEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTEST_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTEST_SRC)))

LIBTEST_DEP := \
		$(LIBTEST_OBJ:.o=.d)

LIBTEST     := $(DIR)/lib$(MODNAME)$(LIBEXT)

TEST_SRC := \
		$(DIR)/test_ckm.cpp \
		$(DIR)/test_cast_model.cpp \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_goldstones.cpp \
		$(DIR)/test_linalg2.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_problems.cpp \
		$(DIR)/test_pv.cpp \
		$(DIR)/test_rk.cpp \
		$(DIR)/test_root_finder.cpp \
		$(DIR)/test_sminput.cpp \
		$(DIR)/test_slha_io.cpp \
		$(DIR)/test_wrappers.cpp

ifneq ($(findstring lattice,$(ALGORITHMS)),)
TEST_SRC +=
endif

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_mssm_solver.cpp \
		$(DIR)/test_two_scale_mssm_initial_guesser.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_loopfunctions.cpp \
		$(DIR)/test_sfermions.cpp \
		$(DIR)/test_CMSSM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSM_higgs_iteration.cpp \
		$(DIR)/test_CMSSM_initial_guesser.cpp \
		$(DIR)/test_CMSSM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_CMSSM_model.cpp \
		$(DIR)/test_CMSSM_spectrum.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-SoftsusyNMSSM --with-CMSSM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_benchmark.cpp \
		$(DIR)/test_CMSSM_slha_output.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-SoftsusyFlavourMSSM --with-MSSM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_MSSM_low_scale_constraint_flavour.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyNMSSM --with-NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_NMSSM_beta_functions.cpp \
		$(DIR)/test_NMSSM_ewsb.cpp \
		$(DIR)/test_NMSSM_high_scale_constraint.cpp \
		$(DIR)/test_NMSSM_initial_guesser.cpp \
		$(DIR)/test_NMSSM_low_scale_constraint.cpp \
		$(DIR)/test_NMSSM_one_loop_spectrum.cpp \
		$(DIR)/test_NMSSM_self_energies.cpp \
		$(DIR)/test_NMSSM_spectrum.cpp \
		$(DIR)/test_NMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_NMSSM_tree_level_spectrum.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyNMSSM --with-SMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_SMSSM_beta_functions.cpp \
		$(DIR)/test_SMSSM_ewsb.cpp \
		$(DIR)/test_SMSSM_one_loop_spectrum.cpp \
		$(DIR)/test_SMSSM_tree_level_spectrum.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-sm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-sm --with-smcw),yes yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm_smcw_integration.cpp
endif
endif
ifeq ($(shell $(FSCONFIG) --with-CMSSM --with-NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_NMSSM_linking.cpp
endif

ifeq ($(shell $(FSCONFIG) --with-CMSSM --with-CMSSMNoFV),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMNoFV_beta_functions.cpp \
		$(DIR)/test_CMSSMNoFV_tree_level_spectrum.cpp \
		$(DIR)/test_CMSSMNoFV_low_scale_constraint.cpp
endif

TEST_SH := \
		$(DIR)/test_space_dir.sh

ifeq ($(ENABLE_LOOPTOOLS),yes)
TEST_SH +=	$(DIR)/test_pv_crosschecks.sh

TEST_PV_EXE := \
		$(DIR)/test_pv_fflite.x \
		$(DIR)/test_pv_looptools.x \
		$(DIR)/test_pv_softsusy.x

$(DIR)/test_pv_crosschecks.sh.log: $(TEST_PV_EXE)

ifneq (,$(findstring test,$(MAKECMDGOALS)))
ALLDEP += $(LIBFFLITE_DEP)
endif
endif

ifeq ($(shell $(FSCONFIG) --with-CMSSM),yes)
TEST_SH += \
		$(DIR)/test_standalone.sh \
		$(DIR)/test_run_examples.sh \
		$(DIR)/test_CMSSM_slha_doubled_blocks.sh
TEST_SRC += \
		$(DIR)/test_CMSSM_slha.cpp \
		$(DIR)/test_CMSSM_slha_input.cpp \
		$(DIR)/test_CMSSM_info.cpp
endif

ifeq ($(shell $(FSCONFIG) --with-SM),yes)
TEST_SRC += \
		$(DIR)/test_SM_beta_functions.cpp \
		$(DIR)/test_SM_low_scale_constraint.cpp
endif

ifeq ($(shell $(FSCONFIG) --with-NSM),yes)
TEST_SRC += \
		$(DIR)/test_NSM_low_scale_constraint.cpp
endif

ifeq ($(shell $(FSCONFIG) --with-lowMSSM --with-CMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_lowMSSM.sh
endif

TEST_META := \
		$(DIR)/test_BetaFunction.m \
		$(DIR)/test_CConversion.m \
		$(DIR)/test_Constraint.m \
		$(DIR)/test_EWSB.m \
		$(DIR)/test_Parameters.m \
		$(DIR)/test_TreeMasses.m \
		$(DIR)/test_SelfEnergies.m \
		$(DIR)/test_TextFormatting.m \
		$(DIR)/test_ThresholdCorrections.m \
		$(DIR)/test_Vertices.m

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_EXE_LOG  := $(TEST_EXE:.x=.x.log)

TEST_SH_LOG   := $(TEST_SH:.sh=.sh.log)

TEST_META_LOG := $(TEST_META:.m=.m.log)

TEST_LOG      := $(TEST_EXE_LOG) $(TEST_SH_LOG)
ifeq ($(ENABLE_META),yes)
TEST_LOG      += $(TEST_META_LOG)
endif

ifeq ($(ENABLE_LOOPTOOLS),yes)
TEST_EXE += $(TEST_PV_EXE)

$(DIR)/test_pv_fflite.x    : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_FFLITE
$(DIR)/test_pv_looptools.x : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_LOOPTOOLS
$(DIR)/test_pv_softsusy.x  : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_SOFTSUSY
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-log \
		execute-tests execute-meta-tests execute-compiled-tests

all-$(MODNAME): $(LIBTEST) $(TEST_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(TEST_DEP)
		-rm -f $(LIBTEST_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(TEST_OBJ)
		-rm -f $(LIBTEST_OBJ)

clean-$(MODNAME)-log:
		-rm -f $(TEST_LOG)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
                  clean-$(MODNAME)-log
		-rm -f $(LIBTEST)
		-rm -f $(TEST_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

$(DIR)/%.x.log: $(DIR)/%.x
		@rm -f $@
		@echo "**************************************************" >> $@;
		@echo "* executing test: $< " >> $@;
		@echo "**************************************************" >> $@;
		@$< --log_level=test_suite >> $@ 2>&1; \
		if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

$(DIR)/%.m.log: $(DIR)/%.m $(META_SRC)
		@rm -f $@
		@echo "**************************************************" >> $@;
		@echo "* executing test: $< " >> $@;
		@echo "**************************************************" >> $@;
		@"$(MATH)" -run "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; \
		Quit[TestSuite\`GetNumberOfFailedTests[]]" >> $@ 2>&1; \
		if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

$(DIR)/%.sh.log: $(DIR)/%.sh
		@rm -f $@
		@echo "**************************************************" >> $@;
		@echo "* executing test: $< " >> $@;
		@echo "**************************************************" >> $@;
		@$< >> $@ 2>&1; \
		if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

execute-tests:  $(TEST_LOG)

ifeq ($(ENABLE_META),yes)
execute-meta-tests: $(TEST_META_LOG)
else
execute-meta-tests:
endif

execute-compiled-tests: $(TEST_EXE_LOG)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_lowMSSM.sh.log: $(RUN_CMSSM_EXE) $(RUN_lowMSSM_EXE)

$(DIR)/test_cast_model.x: $(DIR)/test_cast_model.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_ckm.x: $(DIR)/test_ckm.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_betafunction.x: $(DIR)/test_betafunction.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_goldstones.x: $(DIR)/test_goldstones.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_linalg2.x: $(DIR)/test_linalg2.o
		$(CXX) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(DIR)/test_minimizer.x: $(DIR)/test_minimizer.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_pv.x: $(DIR)/test_pv.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_rk.x: $(DIR)/test_rk.o $(LIBLEGACY) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_root_finder.x: $(DIR)/test_root_finder.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_sminput.x: $(DIR)/test_sminput.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_slha_io.x: $(DIR)/test_slha_io.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(DIR)/test_wrappers.x: $(DIR)/test_wrappers.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_initial_guesser.x: $(DIR)/test_two_scale_mssm_initial_guesser.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_two_scale_sm_smcw_integration.x: $(DIR)/test_two_scale_sm_smcw_integration.o $(LIBsmcw) $(LIBsm) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_two_scale_sm.x: $(DIR)/test_two_scale_sm.o $(LIBsm) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_CMSSM_NMSSM_linking.x: $(DIR)/test_CMSSM_NMSSM_linking.o $(LIBCMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ifeq ($(ENABLE_LOOPTOOLS),yes)
$(DIR)/test_pv_fflite.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(LIBFFLITE)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_pv_looptools.x: $(DIR)/test_pv_crosschecks.cpp $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(LOOPTOOLSLIBS) $(BOOSTTESTLIBS) $(FLIBS)

$(DIR)/test_pv_softsusy.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(filter-out %pv.o,$(LIBFLEXI_OBJ)) $(LIBLEGACY)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)
endif

$(DIR)/test_benchmark.x: $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_loopfunctions.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_sfermions.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_model.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_info.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_initial_guesser.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_higgs_iteration.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_low_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_MSSM_low_scale_constraint_flavour.x: $(LIBSoftsusyFlavourMSSM) $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha_output.x: $(DIR)/test_CMSSM_slha_output.o $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(EXAMPLES_EXE) $(DIR)/test_CMSSM_slha_output.in.spc
		$(CXX) -o $@ $(call abspathx,$< $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_CMSSM_slha_input.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_initial_guesser.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_low_scale_constraint.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_self_energies.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMNoFV_beta_functions.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_CMSSMNoFV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_CMSSMNoFV_low_scale_constraint.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SM_beta_functions.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_low_scale_constraint.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NSM_low_scale_constraint.x: $(LIBNSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

# general test rule which links all libraries needed for a generated model
$(DIR)/test_%.x: $(DIR)/test_%.o
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

# add boost and eigen flags for the test object files and dependencies
$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBTEST): $(LIBTEST_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBTEST): $(LIBTEST_OBJ)
		$(MAKELIB) $@ $^ $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)
endif

ALLDEP += $(LIBTEST_DEP) $(TEST_DEP)
ALLLIB += $(LIBTEST)
ALLTST += $(TEST_EXE)
