DIR      := test
MODNAME  := test
WITH_$(MODNAME) := yes

LIBTEST_SRC := \
		$(DIR)/run_cmd.cpp \
		$(DIR)/stopwatch.cpp

LIBTEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTEST_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTEST_SRC)))

LIBTEST_DEP := \
		$(LIBTEST_OBJ:.o=.d)

LIBTEST     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

TEST_SRC := \
		$(DIR)/test_array_view.cpp \
		$(DIR)/test_ckm.cpp \
		$(DIR)/test_cast_model.cpp \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_lowe.cpp \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_derivative.cpp \
		$(DIR)/test_effective_couplings.cpp \
		$(DIR)/test_eigen_utils.cpp \
		$(DIR)/test_ewsb_solver.cpp \
		$(DIR)/test_fixed_point_iterator.cpp \
		$(DIR)/test_goldstones.cpp \
		$(DIR)/test_gsl_vector.cpp \
		$(DIR)/test_linalg2.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_MSSM_2L_limits.cpp \
		$(DIR)/test_namespace_collisions.cpp \
		$(DIR)/test_numerics.cpp \
		$(DIR)/test_problems.cpp \
		$(DIR)/test_pv.cpp \
		$(DIR)/test_raii.cpp \
		$(DIR)/test_rk.cpp \
		$(DIR)/test_root_finder.cpp \
		$(DIR)/test_scan.cpp \
		$(DIR)/test_sminput.cpp \
		$(DIR)/test_slha_io.cpp \
		$(DIR)/test_sum.cpp \
		$(DIR)/test_thread_pool.cpp \
		$(DIR)/test_threshold_loop_functions.cpp \
		$(DIR)/test_which.cpp \
		$(DIR)/test_wrappers.cpp

TEST_SH := \
		$(DIR)/test_depgen.sh \
		$(DIR)/test_run_examples.sh \
		$(DIR)/test_run_all_spectrum_generators.sh \
		$(DIR)/test_space_dir.sh

ifneq ($(findstring lattice,$(ALGORITHMS)),)
TEST_SRC +=
endif

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp
ifeq ($(WITH_SoftsusyMSSM),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_mssm_solver.cpp \
		$(DIR)/test_two_scale_mssm_initial_guesser.cpp
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_loopfunctions.cpp \
		$(DIR)/test_sfermions.cpp \
		$(DIR)/test_CMSSM_beta_function_benchmark.cpp \
		$(DIR)/test_CMSSM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSM_higgs_iteration.cpp \
		$(DIR)/test_CMSSM_initial_guesser.cpp \
		$(DIR)/test_CMSSM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_CMSSM_model.cpp \
		$(DIR)/test_CMSSM_spectrum.cpp \
		$(DIR)/test_CMSSM_weinberg_angle.cpp
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_CMSSMMassWInput),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMMassWInput_spectrum.cpp
endif

ifeq ($(WITH_CMSSMLowPrecision),yes)
TEST_SRC += \
		$(DIR)/test_CMSSMLowPrecision.cpp
endif
ifeq ($(WITH_SoftsusyMSSM) $(WITH_SoftsusyNMSSM) $(WITH_CMSSM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_benchmark.cpp \
		$(DIR)/test_CMSSM_slha_output.cpp
TEST_SH += \
		$(DIR)/test_CMSSM_gluino.sh
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_SoftsusyFlavourMSSM) $(WITH_CMSSMCKM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCKM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_tree_level_spectrum.cpp
endif

ifeq ($(WITH_SoftsusyNMSSM) $(WITH_NMSSM),yes yes)
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

ifeq ($(WITH_SoftsusyMSSM) $(WITH_SoftsusyNMSSM) $(WITH_NMSSM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_NMSSM_benchmark.cpp \
		$(DIR)/test_NMSSM_slha_output.cpp
endif

ifeq ($(WITH_SoftsusyNMSSM) $(WITH_SMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_SMSSM_beta_functions.cpp \
		$(DIR)/test_SMSSM_ewsb.cpp \
		$(DIR)/test_SMSSM_one_loop_spectrum.cpp \
		$(DIR)/test_SMSSM_tree_level_spectrum.cpp
endif
ifeq ($(WITH_sm),yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm.cpp
endif
ifeq ($(WITH_sm) $(WITH_smcw),yes yes)
TEST_SRC += \
		$(DIR)/test_two_scale_sm_smcw_integration.cpp
endif
endif
ifeq ($(WITH_CMSSM) $(WITH_NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_NMSSM_linking.cpp
endif

ifeq ($(WITH_CMSSMNoFV),yes)
TEST_SRC += \
		$(DIR)/test_CMSSMNoFV_two_loop_spectrum.cpp
endif

ifeq ($(WITH_GM2Calc),yes)
TEST_SRC += \
		$(DIR)/test_gm2calc.cpp \
		$(DIR)/test_MSSMNoFV_onshell.cpp
endif

ifeq ($(WITH_GM2Calc) $(WITH_CMSSMNoFV),yes yes)
TEST_SH += \
		$(DIR)/test_CMSSMNoFV_GM2Calc.sh
endif

ifeq ($(WITH_CMSSM) $(WITH_CMSSMNoFV),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMNoFV_beta_functions.cpp \
		$(DIR)/test_CMSSMNoFV_tree_level_spectrum.cpp \
		$(DIR)/test_CMSSMNoFV_low_scale_constraint.cpp
endif

ifeq ($(WITH_CMSSM) $(WITH_cCMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_cCMSSM.sh
endif

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

ifeq ($(WITH_BLSMlightZp),yes)
TEST_SH += \
		$(DIR)/test_BLSMlightZp_ZZp_mixing.sh
endif

ifeq ($(WITH_MSSM),yes)
TEST_SH += \
		$(DIR)/test_MSSM_stable_ewsb_failure.sh \
		$(DIR)/test_standalone.sh
endif

ifeq ($(WITH_CMSSM),yes)
TEST_SH += \
		$(DIR)/test_CMSSM_slha_doubled_blocks.sh \
		$(DIR)/test_CMSSM_memory_leaks.sh \
		$(DIR)/test_CMSSM_profile.sh \
		$(DIR)/test_CMSSM_QedQcd_exception.sh \
		$(DIR)/test_CMSSM_QedQcd_no_convergence.sh
TEST_SRC += \
		$(DIR)/test_CMSSM_database.cpp \
		$(DIR)/test_CMSSM_slha.cpp \
		$(DIR)/test_CMSSM_slha_input.cpp \
		$(DIR)/test_CMSSM_two_loop_spectrum.cpp \
		$(DIR)/test_CMSSM_info.cpp
endif

ifeq ($(WITH_NMSSM),yes)
TEST_SH += \
		$(DIR)/test_NMSSM_profile.sh
endif

ifeq ($(WITH_SM),yes)
TEST_SRC += \
		$(DIR)/test_SM_beta_functions.cpp \
		$(DIR)/test_SM_effective_couplings.cpp \
		$(DIR)/test_SM_gmm2.cpp \
		$(DIR)/test_SM_low_scale_constraint.cpp \
		$(DIR)/test_SM_one_loop_spectrum.cpp \
		$(DIR)/test_SM_higgs_loop_corrections.cpp \
		$(DIR)/test_SM_tree_level_spectrum.cpp \
		$(DIR)/test_SM_three_loop_spectrum.cpp \
		$(DIR)/test_SM_two_loop_spectrum.cpp \
		$(DIR)/test_SM_weinberg_angle.cpp
endif

ifeq ($(WITH_SMHighPrecision),yes)
TEST_SRC += \
		$(DIR)/test_SMHighPrecision_two_loop_spectrum.cpp
endif

ifeq ($(WITH_NSM),yes)
TEST_SRC += \
		$(DIR)/test_NSM_low_scale_constraint.cpp
endif

ifeq ($(WITH_lowMSSM) $(WITH_CMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_lowMSSM.sh
endif

ifeq ($(WITH_lowNMSSM) $(WITH_SoftsusyMSSM) $(WITH_SoftsusyNMSSM),yes yes yes)
TEST_SH += \
		$(DIR)/test_lowNMSSM_spectrum.sh
endif

ifeq ($(WITH_CMSSMGSLHybrid) $(WITH_CMSSMGSLHybridS) $(WITH_CMSSMGSLBroyden) $(WITH_CMSSMGSLNewton) $(WITH_CMSSMFPIRelative) $(WITH_CMSSMFPIAbsolute) $(WITH_CMSSMFPITadpole),yes yes yes yes yes yes yes)
TEST_SRC += \
		$(DIR)/test_compare_ewsb_solvers.cpp
endif

ifeq ($(WITH_CMSSMCPV),yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCPV_ewsb.cpp
endif

ifeq ($(WITH_CMSSMCPV) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCPV_tree_level_spectrum.cpp
TEST_SH += \
		$(DIR)/test_CMSSMCPV_spectrum.sh
endif

ifeq ($(WITH_NMSSMCPV),yes)
TEST_SRC += \
		$(DIR)/test_NMSSMCPV_ewsb.cpp
endif

ifeq ($(WITH_NMSSMCPV) $(WITH_NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_NMSSMCPV_tree_level_spectrum.cpp
endif

ifeq ($(WITH_NUTNMSSM),yes)
TEST_SH += \
		$(DIR)/test_NUTNMSSM.sh
endif

ifeq ($(WITH_NUTNMSSM) $(WITH_SoftsusyNMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_NUTNMSSM_spectrum.cpp
endif

ifeq ($(WITH_HSSUSY),yes)
TEST_SH += \
		$(DIR)/test_HSSUSY_SUSYHD.sh
endif

ifeq ($(WITH_NUHMSSMalttower) $(WITH_NUHMSSMalt),yes yes)
TEST_SH += \
		$(DIR)/test_NUHMSSMalttower.sh
endif

ifeq ($(WITH_HSSUSY) $(WITH_MSSMtower) $(WITH_MSSMMuBMu),yes yes yes)
TEST_SH += \
		$(DIR)/test_MSSMtower.sh
endif

ifeq ($(WITH_MSSMtower),yes)
TEST_SH += \
		$(DIR)/test_MSSMtower_librarylink.sh
endif

ifeq ($(WITH_MSSMtower) $(WITH_MSSMNoFVtower),yes yes)
TEST_SH += \
		$(DIR)/test_MSSMNoFVtower.sh
endif

ifeq ($(WITH_NMSSMtower) $(WITH_lowNMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_NMSSMtower.sh
endif

ifeq ($(WITH_SMHighPrecision) $(WITH_SMtower),yes yes)
TEST_SH += \
		$(DIR)/test_SMtower.sh
endif

ifeq ($(WITH_SM) $(WITH_SMtower),yes yes)
TEST_META += \
		$(DIR)/test_multiple_librarylinks.m
endif

ifeq ($(WITH_SM),yes)
TEST_SH += \
		$(DIR)/test_flexiblesusy-config.sh
endif

ifeq ($(WITH_THDMIIMSSMBC) $(WITH_THDMIIMSSMBCApprox) $(WITH_HGTHDMIIMSSMBC) $(WITH_HGTHDMIIMSSMBCApprox),yes yes yes yes)

TEST_SH += \
		$(DIR)/test_THDMIIMSSMBCFull_approximation.sh
endif


ifeq ($(shell $(FSCONFIG) --with-HGTHDMIIMSSMBC --with-MSSM),yes yes)
TEST_META += \
		$(DIR)/test_HGTHDM_threshold_corrections_scale_invariance.m
endif

ifeq ($(shell $(FSCONFIG) --with-THDMIIMSSMBC --with-MSSM),yes yes)
TEST_META += \
		$(DIR)/test_THDM_threshold_corrections_scale_invariance.m
endif

TEST_META := \
		$(DIR)/test_BetaFunction.m \
		$(DIR)/test_CConversion.m \
		$(DIR)/test_Constraint.m \
		$(DIR)/test_EWSB.m \
		$(DIR)/test_HSSUSY_thresholds.m \
		$(DIR)/test_LoopFunctions.m \
		$(DIR)/test_MSSM_2L_analytic.m \
		$(DIR)/test_Parameters.m \
		$(DIR)/test_ReadSLHA.m \
		$(DIR)/test_RGIntegrator.m \
		$(DIR)/test_SelfEnergies.m \
		$(DIR)/test_TextFormatting.m \
		$(DIR)/test_THDM_threshold_corrections.m \
		$(DIR)/test_THDM_threshold_corrections_gauge.m \
		$(DIR)/test_ThreeLoopQCD.m \
		$(DIR)/test_ThresholdCorrections.m \
		$(DIR)/test_TreeMasses.m \
		$(DIR)/test_Vertices.m

ifeq ($(WITH_SM),yes)
TEST_META += \
		$(DIR)/test_SM_3loop_beta.m
endif

ifeq ($(WITH_SMnoGUT),yes)
TEST_META += \
		$(DIR)/test_SMnoGUT_3loop_beta.m
endif

ifeq ($(WITH_CMSSM),yes)
TEST_META += \
		$(DIR)/test_CMSSM_3loop_beta.m \
		$(DIR)/test_CMSSM_librarylink.m \
		$(DIR)/test_CMSSM_librarylink_parallel.m
TEST_SH += \
		$(DIR)/test_CMSSM_librarylink_slha.sh
endif

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
		clean-$(MODNAME)-dep clean-$(MODNAME)-log \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj \
		execute-tests execute-meta-tests execute-compiled-tests \
		execute-shell-tests

all-$(MODNAME): $(LIBTEST) $(TEST_EXE) $(TEST_LOG)
		@true

clean-$(MODNAME)-dep:
		-rm -f $(TEST_DEP)
		-rm -f $(LIBTEST_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBTEST)

clean-$(MODNAME)-obj:
		-rm -f $(TEST_OBJ)
		-rm -f $(LIBTEST_OBJ)

clean-$(MODNAME)-log:
		-rm -f $(TEST_LOG)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
                  clean-$(MODNAME)-lib clean-$(MODNAME)-log
		-rm -f $(TEST_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

execute-tests:  $(TEST_LOG)

ifeq ($(ENABLE_META),yes)
execute-meta-tests: $(TEST_META_LOG)
else
execute-meta-tests:
endif

execute-compiled-tests: $(TEST_EXE_LOG)

execute-shell-tests: $(TEST_SH_LOG)

PTR = print_test_result() { \
	if [ $$1 = 0 ]; then \
		printf "%-66s %4s\n" "$$2" "OK"; \
	else \
		printf "%-66s %4s\n" "$$2" "FAILED"; \
	fi \
}

$(DIR)/%.x.log: $(DIR)/%.x
		@rm -f $@
		@$(PTR); \
		BOOST_TEST_CATCH_SYSTEM_ERRORS="no" \
		$< --log_level=test_suite >> $@ 2>&1; \
		print_test_result $$? $<

$(DIR)/%.m.log: $(DIR)/%.m $(META_SRC)
		@rm -f $@
		@$(PTR); \
		"$(MATH)" -run "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; \
		Quit[TestSuite\`GetNumberOfFailedTests[]]" >> $@ 2>&1; \
		print_test_result $$? $<

$(DIR)/%.sh.log: $(DIR)/%.sh
		@rm -f $@
		@$(PTR); \
		MATH_CMD="$(MATH)" $< >> $@ 2>&1; \
		print_test_result $$? $<

$(DIR)/test_lowMSSM.sh.log: $(RUN_CMSSM_EXE) $(RUN_lowMSSM_EXE)

$(DIR)/test_cast_model.x: $(DIR)/test_cast_model.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_ckm.x: $(DIR)/test_ckm.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_lowe.x: $(DIR)/test_lowe.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_betafunction.x: $(DIR)/test_betafunction.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_effective_couplings.x: $(DIR)/test_effective_couplings.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_ewsb_solver.x: $(DIR)/test_ewsb_solver.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_fixed_point_iterator.x: $(DIR)/test_fixed_point_iterator.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_goldstones.x: $(DIR)/test_goldstones.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_gsl_vector.x: $(DIR)/test_gsl_vector.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_linalg2.x: $(DIR)/test_linalg2.o
		$(CXX) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(DIR)/test_MSSM_2L_limits.x: $(DIR)/test_MSSM_2L_limits.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_minimizer.x: $(DIR)/test_minimizer.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_numerics.x: $(DIR)/test_numerics.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_pv.x: $(DIR)/test_pv.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_rk.x: $(DIR)/test_rk.o $(LIBLEGACY) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_root_finder.x: $(DIR)/test_root_finder.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_sminput.x: $(DIR)/test_sminput.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_slha_io.x: $(DIR)/test_slha_io.o $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_thread_pool.x: $(DIR)/test_thread_pool.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_threshold_loop_functions.x: $(DIR)/test_threshold_loop_functions.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_wrappers.x: $(DIR)/test_wrappers.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_sum.x: $(DIR)/test_sum.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS)

$(DIR)/test_two_scale_mssm_initial_guesser.x: $(DIR)/test_two_scale_mssm_initial_guesser.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_two_scale_sm_smcw_integration.x: $(DIR)/test_two_scale_sm_smcw_integration.o $(LIBsmcw) $(LIBsm) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_two_scale_sm.x: $(DIR)/test_two_scale_sm.o $(LIBsm) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_CMSSM_NMSSM_linking.x: $(DIR)/test_CMSSM_NMSSM_linking.o $(LIBCMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ifeq ($(ENABLE_LOOPTOOLS),yes)
$(DIR)/test_pv_fflite.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(LIBFFLITE)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_looptools.x: $(DIR)/test_pv_crosschecks.cpp $(LIBFLEXI) $(LIBLEGACY)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(LOOPTOOLSLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_softsusy.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(filter-out %pv.o,$(LIBFLEXI_OBJ)) $(LIBLEGACY)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)
endif

$(DIR)/test_CMSSM_benchmark.x: CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)
$(DIR)/test_CMSSM_benchmark.x: $(DIR)/test_CMSSM_benchmark.cpp $(RUN_CMSSM_EXE) $(RUN_SOFTPOINT_EXE) $(LIBTEST)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$<) $(LIBTEST) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_compare_ewsb_solvers.x: $(LIBCMSSMGSLHybrid) $(LIBCMSSMGSLHybridS) $(LIBCMSSMGSLBroyden) $(LIBCMSSMGSLNewton) $(LIBCMSSMFPIRelative) $(LIBCMSSMFPIAbsolute) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_loopfunctions.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_sfermions.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_database.x: $(DIR)/test_CMSSM_database.o $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(THREADLIBS)

$(DIR)/test_CMSSM_model.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_info.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_two_loop_spectrum.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_beta_function_benchmark.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)

$(DIR)/test_CMSSM_initial_guesser.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_higgs_iteration.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_low_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha_output.x: $(DIR)/test_CMSSM_slha_output.o $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(EXAMPLES_EXE) $(DIR)/test_CMSSM_slha_output.in.spc
		$(CXX) -o $@ $(call abspathx,$< $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_CMSSM_slha_input.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCKM_high_scale_constraint.x \
$(DIR)/test_CMSSMCKM_low_scale_constraint.x \
$(DIR)/test_CMSSMCKM_tree_level_spectrum.x: \
	CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

$(DIR)/test_CMSSMCKM_high_scale_constraint.x \
$(DIR)/test_CMSSMCKM_low_scale_constraint.x \
$(DIR)/test_CMSSMCKM_tree_level_spectrum.x: \
	$(LIBSoftsusyFlavourMSSM) $(LIBSoftsusyMSSM) $(LIBCMSSMCKM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_weinberg_angle.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMMassWInput_spectrum.x: $(LIBSoftsusyMSSM) $(LIBCMSSMMassWInput) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMLowPrecision.x: $(LIBSoftsusyMSSM) $(LIBCMSSMLowPrecision) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPV_ewsb.x: $(LIBCMSSMCPV) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMCPV) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSMCPV_ewsb.x: $(LIBNMSSMCPV) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSMCPV_tree_level_spectrum.x: $(LIBNMSSM) $(LIBNMSSMCPV) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)

$(DIR)/test_NMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_initial_guesser.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_low_scale_constraint.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_self_energies.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_benchmark.x: CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)
$(DIR)/test_NMSSM_benchmark.x: $(DIR)/test_NMSSM_benchmark.cpp $(RUN_NMSSM_EXE) $(RUN_SOFTPOINT_EXE) $(LIBTEST)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$<) $(LIBTEST) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_NMSSM_slha_output.x: $(DIR)/test_NMSSM_slha_output.o $(LIBNMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(EXAMPLES_EXE) $(DIR)/test_NMSSM_slha_output.in.spc
		$(CXX) -o $@ $(call abspathx,$< $(LIBNMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_SMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NUTNMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNUTNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMNoFV_beta_functions.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_CMSSMNoFV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_CMSSMNoFV_two_loop_spectrum.x: $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_CMSSMNoFV_low_scale_constraint.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_gm2calc.x: $(LIBMSSMNoFVSLHA2) $(LIBGM2Calc) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_MSSMNoFV_onshell.x: $(LIBMSSMNoFVSLHA2) $(LIBGM2Calc) $(LIBFLEXI) $(LIBLEGACY)

$(DIR)/test_SM_beta_functions.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_effective_couplings.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_gmm2.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_higgs_loop_corrections.x: $(DIR)/test_SM_higgs_loop_corrections.o $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$< $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_SM_low_scale_constraint.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_tree_level_spectrum.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_one_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_three_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_two_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_weinberg_angle.x: $(LIBSoftsusyMSSM) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMHighPrecision_two_loop_spectrum.x: $(LIBSMHighPrecision) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NSM_low_scale_constraint.x: $(LIBNSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

# general test rule which links all libraries needed for a generated model
$(DIR)/test_%.x: $(DIR)/test_%.o
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

# add boost and eigen flags for the test object files and dependencies
$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIBTEST): $(LIBTEST_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)
else
$(LIBTEST): $(LIBTEST_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIBTEST_DEP) $(TEST_DEP)
ALLLIB += $(LIBTEST)
ALLTST += $(TEST_EXE)
