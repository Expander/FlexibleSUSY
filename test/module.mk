DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_linalg2.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_problems.cpp \
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
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-MSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_loopfunctions.cpp \
		$(DIR)/test_MSSM_high_scale_constraint.cpp \
		$(DIR)/test_MSSM_higgs_iteration.cpp \
		$(DIR)/test_MSSM_initial_guesser.cpp \
		$(DIR)/test_MSSM_low_scale_constraint.cpp \
		$(DIR)/test_MSSM_susy_scale_constraint.cpp \
		$(DIR)/test_MSSM_model.cpp \
		$(DIR)/test_MSSM_spectrum.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyMSSM --with-SoftsusyNMSSM --with-MSSM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_benchmark.cpp \
		$(DIR)/test_MSSM_slha_output.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-SoftsusyNMSSM --with-NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_NMSSM_beta_functions.cpp \
		$(DIR)/test_NMSSM_ewsb.cpp \
		$(DIR)/test_NMSSM_high_scale_constraint.cpp \
		$(DIR)/test_NMSSM_initial_guesser.cpp \
		$(DIR)/test_NMSSM_low_scale_constraint.cpp \
		$(DIR)/test_NMSSM_one_loop_spectrum.cpp \
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
ifeq ($(shell $(FSCONFIG) --with-MSSM --with-NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_MSSM_NMSSM_linking.cpp
endif

TEST_SH := \
		test/test_space_dir.sh

ifeq ($(shell $(FSCONFIG) --with-lowMSSM --with-MSSM),yes yes)
TEST_SH += \
		test/test_lowMSSM.sh
endif

TEST_META := \
		test/test_BetaFunction.m \
		test/test_CConversion.m \
		test/test_Constraint.m \
		test/test_EWSB.m \
		test/test_Parameters.m \
		test/test_TreeMasses.m \
		test/test_SelfEnergies.m \
		test/test_TextFormatting.m \
		test/test_ThresholdCorrections.m

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_EXE_LOG  := $(TEST_EXE:.x=.x.log)

TEST_SH_LOG   := $(TEST_SH:.sh=.sh.log)

TEST_META_LOG := $(TEST_META:.m=.m.log)

TEST_LOG      := $(TEST_EXE_LOG) $(TEST_SH_LOG) $(TEST_META_LOG)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-log \
		execute-tests execute-meta-tests execute-compiled-tests

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME)-log:
		-rm -f $(TEST_LOG)

clean-$(MODNAME):
		-rm -f $(TEST_OBJ)
		-rm -f $(TEST_LOG)

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(TEST_DEP)
		-rm -f $(TEST_EXE)

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
		@$(MATH) -run "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; \
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

execute-meta-tests: $(TEST_META_LOG)

execute-compiled-tests: $(TEST_EXE_LOG)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_lowMSSM.sh.log: $(RUN_MSSM_EXE) $(RUN_lowMSSM_EXE)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_betafunction.x: $(DIR)/test_betafunction.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_linalg2.x: $(DIR)/test_linalg2.o
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(LAPACKLIBS) $(FLIBS)

$(DIR)/test_minimizer.x: $(DIR)/test_minimizer.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_rk.x: $(DIR)/test_rk.o $(LIBLEGACY) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_root_finder.x: $(DIR)/test_root_finder.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_sminput.x: $(DIR)/test_sminput.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_slha_io.x: $(DIR)/test_slha_io.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_wrappers.x: $(DIR)/test_wrappers.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_initial_guesser.x: $(DIR)/test_two_scale_mssm_initial_guesser.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm_smcw_integration.x: $(DIR)/test_two_scale_sm_smcw_integration.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm.x: $(DIR)/test_two_scale_sm.o $(LIBSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_MSSM_NMSSM_linking.x: $(DIR)/test_MSSM_NMSSM_linking.o $(LIBMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_benchmark.x: $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_loopfunctions.x: $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_model.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_initial_guesser.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_higgs_iteration.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_low_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSM_slha_output.x: $(DIR)/test_MSSM_slha_output.o $(LIBMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(EXAMPLES_EXE) $(DIR)/test_MSSM_slha_output.in.spc
		$(CXX) -o $@ $< $(LIBMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_MSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_initial_guesser.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_low_scale_constraint.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_NMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)

# general test rule which links all libraries needed for a generated model
$(DIR)/test_%.x: $(DIR)/test_%.o
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

# add boost and eigen flags for the test object files and dependencies
$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
