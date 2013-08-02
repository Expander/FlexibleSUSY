DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_rk.cpp \
		$(DIR)/test_root_finder.cpp \
		$(DIR)/test_wrappers.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp
ifeq ($(shell $(FSCONFIG) --with-smssm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_mssm_solver.cpp \
		$(DIR)/test_two_scale_mssm_initial_guesser.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-smssm --with-MSSM),yes yes)
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
ifeq ($(shell $(FSCONFIG) --with-sm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm.cpp
endif
ifeq ($(shell $(FSCONFIG) --with-sm --with-smcw),yes yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm_smcw_integration.cpp
endif
endif

TEST_META := \
		test/test_CConversion.m \
		test/test_Constraint.m \
		test/test_TreeMasses.m \
		test/test_SelfEnergies.m \
		test/test_TextFormatting.m

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_SCRIPT := \
		$(DIR)/execute_test.sh

TEST_META_SCRIPT := \
		$(DIR)/execute_meta_code_test.sh

TEST_LOG := \
		$(TEST_EXE:.x=.x.log) \
		$(TEST_META:.m=.m.log)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME)-log:
		rm -rf $(TEST_LOG)

clean-$(MODNAME):
		rm -rf $(TEST_OBJ)
		rm -rf $(TEST_LOG)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(TEST_DEP)
		rm -rf $(TEST_EXE)

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

execute-tests:  $(TEST_LOG)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_betafunction.d $(DIR)/test_betafunction.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_betafunction.x: $(DIR)/test_betafunction.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_minimizer.d $(DIR)/test_minimizer.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_minimizer.x: $(DIR)/test_minimizer.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_rk.d $(DIR)/test_rk.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_rk.x: $(DIR)/test_rk.o $(LIBLEGACY) $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_root_finder.d $(DIR)/test_root_finder.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_root_finder.x: $(DIR)/test_root_finder.o $(LIBFLEXI)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS)

$(DIR)/test_wrappers.d $(DIR)/test_wrappers.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_wrappers.x: $(DIR)/test_wrappers.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_mssm_initial_guesser.d $(DIR)/test_two_scale_mssm_initial_guesser.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_two_scale_mssm_initial_guesser.x: $(DIR)/test_two_scale_mssm_initial_guesser.o $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(FLIBS) $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm_smcw_integration.d $(DIR)/test_two_scale_sm_smcw_integration.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_two_scale_sm_smcw_integration.x: $(DIR)/test_two_scale_sm_smcw_integration.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_sm.x: $(DIR)/test_two_scale_sm.o $(LIBSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

$(DIR)/test_loopfunctions.d $(DIR)/test_loopfunctions.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_loopfunctions.x: $(DIR)/test_loopfunctions.o $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_model.d $(DIR)/test_MSSM_model.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_model.x: $(DIR)/test_MSSM_model.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_initial_guesser.d $(DIR)/test_MSSM_initial_guesser.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_initial_guesser.x: $(DIR)/test_MSSM_initial_guesser.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_higgs_iteration.d $(DIR)/test_MSSM_higgs_iteration.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_higgs_iteration.x: $(DIR)/test_MSSM_higgs_iteration.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_high_scale_constraint.d $(DIR)/test_MSSM_high_scale_constraint.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_high_scale_constraint.x: $(DIR)/test_MSSM_high_scale_constraint.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_low_scale_constraint.d $(DIR)/test_MSSM_low_scale_constraint.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_low_scale_constraint.x: $(DIR)/test_MSSM_low_scale_constraint.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_susy_scale_constraint.d $(DIR)/test_MSSM_susy_scale_constraint.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_susy_scale_constraint.x: $(DIR)/test_MSSM_susy_scale_constraint.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

$(DIR)/test_MSSM_spectrum.d $(DIR)/test_MSSM_spectrum.x: CPPFLAGS += $(EIGENFLAGS)
$(DIR)/test_MSSM_spectrum.x: $(DIR)/test_MSSM_spectrum.o $(LIBSMSSM) $(LIBMSSM) $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS) $(GSLLIBS) $(FLIBS) $(LOOPTOOLSLIBS)

%.x: %.o $(ALLLIB)
		$(CXX) -o $@ $^ $(BOOSTTESTLIBS)

# add boost flags for the test object files only
$(TEST_OBJ): CPPFLAGS += $(BOOSTFLAGS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
