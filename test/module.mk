DIR      := test
MODNAME  := test
WITH_$(MODNAME) := yes

LIBTEST_SRC := \
		$(DIR)/error_count.cpp \
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
		$(DIR)/test_cast_model.cpp \
		$(DIR)/test_ckm.cpp \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_derivative.cpp \
		$(DIR)/test_effective_couplings.cpp \
		$(DIR)/test_eigen_utils.cpp \
		$(DIR)/test_ewsb_solver.cpp \
		$(DIR)/test_fixed_point_iterator.cpp \
		$(DIR)/test_goldstones.cpp \
		$(DIR)/test_gsl_vector.cpp \
		$(DIR)/test_linalg2.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_mssm_twoloop_mb.cpp \
		$(DIR)/test_mssm_twoloop_mt.cpp \
		$(DIR)/test_MSSM_2L_limits.cpp \
		$(DIR)/test_namespace_collisions.cpp \
		$(DIR)/test_numerics.cpp \
		$(DIR)/test_problems.cpp \
		$(DIR)/test_pv.cpp \
		$(DIR)/test_raii.cpp \
		$(DIR)/test_root_finder.cpp \
		$(DIR)/test_scan.cpp \
		$(DIR)/test_sminput.cpp \
		$(DIR)/test_slha_io.cpp \
		$(DIR)/test_sum.cpp \
		$(DIR)/test_thread_pool.cpp \
		$(DIR)/test_threshold_corrections.cpp \
		$(DIR)/test_threshold_loop_functions.cpp \
		$(DIR)/test_spectrum_generator_settings.cpp \
		$(DIR)/test_which.cpp \
		$(DIR)/test_wrappers.cpp

TEST_SH := \
		$(DIR)/test_depgen.sh \
		$(DIR)/test_run_examples.sh \
		$(DIR)/test_run_all_spectrum_generators.sh \
		$(DIR)/test_space_dir.sh

TEST_META := \
		$(DIR)/test_BetaFunction.m \
		$(DIR)/test_CConversion.m \
		$(DIR)/test_Constraint.m \
		$(DIR)/test_EWSB.m \
		$(DIR)/test_LoopFunctions.m \
		$(DIR)/test_MSSM_2L_analytic.m \
		$(DIR)/test_MSSM_2L_mt.m \
		$(DIR)/test_MSSM_2L_yb_softsusy.m \
		$(DIR)/test_MSSM_2L_yt.m \
		$(DIR)/test_MSSM_2L_yt_loopfunction.m \
		$(DIR)/test_MSSM_2L_yt_softsusy.m \
		$(DIR)/test_Parameters.m \
		$(DIR)/test_ReadSLHA.m \
		$(DIR)/test_RGIntegrator.m \
		$(DIR)/test_SelfEnergies.m \
		$(DIR)/test_SemiAnalytic.m \
		$(DIR)/test_TextFormatting.m \
		$(DIR)/test_THDM_threshold_corrections.m \
		$(DIR)/test_THDM_threshold_corrections_gauge.m \
		$(DIR)/test_ThreeLoopQCD.m \
		$(DIR)/test_ThresholdCorrections.m \
		$(DIR)/test_TreeMasses.m \
		$(DIR)/test_Vertices_SortCp.m \
		$(DIR)/test_Vertices_colorsum.m

ifneq ($(findstring lattice,$(SOLVERS)),)
TEST_SRC +=
endif # ifneq($(findstring lattice,$(SOLVERS)),)

ifneq ($(findstring two_scale,$(SOLVERS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp

ifeq ($(WITH_SoftsusyMSSM),yes)
TEST_SRC += \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_legacy_diagonalization.cpp \
		$(DIR)/test_lowe.cpp \
		$(DIR)/test_QedQcd.cpp \
		$(DIR)/test_rk.cpp \
		$(DIR)/test_two_scale_mssm_solver.cpp \
		$(DIR)/test_two_scale_mssm_initial_guesser.cpp
endif

ifeq ($(WITH_SM) $(WITH_SoftsusyMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_SM_weinberg_angle.cpp \
		$(DIR)/test_SM_weinberg_angle_meta.cpp
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_loopfunctions.cpp \
		$(DIR)/test_sfermions.cpp \
		$(DIR)/test_CMSSM_beta_function_benchmark.cpp \
		$(DIR)/test_CMSSM_database.cpp \
		$(DIR)/test_CMSSM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSM_higgs_iteration.cpp \
		$(DIR)/test_CMSSM_initial_guesser.cpp \
		$(DIR)/test_CMSSM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSM_model.cpp \
		$(DIR)/test_CMSSM_spectrum.cpp \
		$(DIR)/test_CMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_CMSSM_weinberg_angle.cpp \
		$(DIR)/test_CMSSM_weinberg_angle_meta.cpp
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
		$(DIR)/test_CMSSM_slha_output.cpp
TEST_SH += \
		$(DIR)/test_CMSSM_gluino.sh
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_SoftsusyFlavourMSSM) $(WITH_CMSSMCKM),yes yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCKM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_tree_level_spectrum.cpp
TEST_SH += \
		$(DIR)/test_CMSSMCKM_spectrum.sh
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

ifeq ($(ENABLE_SQLITE) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_database.cpp
endif

ifeq ($(WITH_MRSSM2),yes)
TEST_SRC += \
		$(DIR)/test_MRSSM2_gmm2.cpp
endif

endif # ifneq ($(findstring two_scale,$(SOLVERS)),)

ifneq ($(findstring semi_analytic,$(SOLVERS)),)

ifeq ($(WITH_CMSSMSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_CMSSMSemiAnalytic_ewsb.cpp \
		$(DIR)/test_CMSSMSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_CMSSMCPVSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCPVSemiAnalytic_ewsb.cpp \
		$(DIR)/test_CMSSMCPVSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_CNMSSM), yes)
TEST_SRC += \
		$(DIR)/test_CNMSSM_ewsb.cpp \
		$(DIR)/test_CNMSSM_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_CE6SSM), yes)
TEST_SRC += \
		$(DIR)/test_CE6SSM_ewsb.cpp \
		$(DIR)/test_CE6SSM_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_lowNUHMSSMSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_lowNUHMSSMSemiAnalytic_ewsb.cpp \
		$(DIR)/test_lowNUHMSSMSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_munuSSMSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_munuSSMSemiAnalytic_ewsb.cpp \
		$(DIR)/test_munuSSMSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_SMSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_SMSemiAnalytic_ewsb.cpp \
		$(DIR)/test_SMSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_SSMSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_SSMSemiAnalytic_ewsb.cpp \
		$(DIR)/test_SSMSemiAnalytic_semi_analytic_solutions.cpp
endif

ifeq ($(WITH_THDMIIEWSBAtMZSemiAnalytic), yes)
TEST_SRC += \
		$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_ewsb.cpp \
		$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions.cpp
endif

ifneq ($(findstring two_scale,$(SOLVERS)),)
ifeq ($(WITH_CMSSM) $(WITH_CMSSMSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_CMSSMSemiAnalytic_spectrum.sh

TEST_SRC += \
		$(DIR)/test_CMSSMSemiAnalytic_consistent_solutions.cpp \
		$(DIR)/test_CMSSMSemiAnalytic_ewsb_solution.cpp
endif

ifeq ($(WITH_CMSSMCPV) $(WITH_CMSSMCPVSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_CMSSMCPVSemiAnalytic_spectrum.sh
endif

ifeq ($(WITH_NMSSM) $(WITH_CNMSSM), yes yes)
TEST_SH += \
		$(DIR)/test_CNMSSM_spectrum.sh

TEST_SRC += \
		$(DIR)/test_CNMSSM_consistent_solutions.cpp
endif

ifeq ($(WITH_E6SSM) $(WITH_CE6SSM), yes yes)
TEST_SH += \
		$(DIR)/test_CE6SSM_spectrum.sh

TEST_SRC += \
		$(DIR)/test_CE6SSM_consistent_solutions.cpp
endif

ifeq ($(WITH_lowNUHMSSM) $(WITH_lowNUHMSSMSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_lowNUHMSSMSemiAnalytic_spectrum.sh

TEST_SRC += \
		$(DIR)/test_lowNUHMSSMSemiAnalytic_consistent_solutions.cpp
endif

ifeq ($(WITH_munuSSM) $(WITH_munuSSMSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_munuSSMSemiAnalytic_spectrum.sh
endif

ifeq ($(WITH_SM) $(WITH_SMSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_SMSemiAnalytic_spectrum.sh

TEST_SRC += \
		$(DIR)/test_SMSemiAnalytic_consistent_solutions.cpp
endif

ifeq ($(WITH_SSM) $(WITH_SSMSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_SSMSemiAnalytic_spectrum.sh

TEST_SRC += \
		$(DIR)/test_SSMSemiAnalytic_consistent_solutions.cpp
endif

ifeq ($(WITH_THDMII) $(WITH_THDMIIEWSBAtMZSemiAnalytic), yes yes)
TEST_SH += \
		$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_spectrum.sh

TEST_SRC += \
		$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_consistent_solutions.cpp
endif

endif # ifneq ($(findstring two_scale,$(SOLVERS)),)

endif # ifneq ($(findstring semi_analytic,$(SOLVERS)),)

ifeq ($(WITH_CMSSM) $(WITH_NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_NMSSM_linking.cpp
endif

ifeq ($(WITH_CMSSM),yes)
TEST_SRC += \
		$(DIR)/test_CMSSM_effective_couplings.cpp
endif

ifeq ($(WITH_CMSSMNoFV),yes)
TEST_SH += \
		$(DIR)/test_CMSSMNoFV_profile.sh
endif

ifeq ($(WITH_SoftsusyMSSM) $(WITH_SoftsusyNMSSM) $(WITH_CMSSMNoFV),yes yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMNoFV_benchmark.cpp \
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

ifeq ($(WITH_SoftsusyMSSM) $(WITH_CMSSM) $(WITH_CMSSMNoFV),yes yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMNoFV_beta_functions.cpp \
		$(DIR)/test_CMSSMNoFV_tree_level_spectrum.cpp \
		$(DIR)/test_CMSSMNoFV_low_scale_constraint.cpp \
		$(DIR)/test_CMSSMNoFV_weinberg_angle_meta.cpp
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
endif # ifeq ($(ENABLE_LOOPTOOLS),yes)

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
		$(DIR)/test_CMSSM_QedQcd_no_convergence.sh \
		$(DIR)/test_CMSSM_streams.sh
TEST_SRC += \
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
		$(DIR)/test_SM_two_loop_spectrum.cpp
endif

ifeq ($(WITH_SMHighPrecision),yes)
TEST_SRC += \
		$(DIR)/test_SMHighPrecision_two_loop_spectrum.cpp
endif

ifeq ($(WITH_SMSU3),yes)
TEST_SRC += \
		$(DIR)/test_SMSU3_low_scale_constraint.cpp
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
TEST_META += \
		$(DIR)/test_CMSSMCPV_librarylink.m
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

ifeq ($(WITH_VCMSSM) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_VCMSSM_ewsb.cpp
TEST_SH += \
		$(DIR)/test_VCMSSM_spectrum.sh
endif

ifeq ($(WITH_HSSUSY),yes)
TEST_SH += \
		$(DIR)/test_HSSUSY_SUSYHD.sh
endif

ifeq ($(WITH_NUHMSSMaltEFTHiggs) $(WITH_NUHMSSMalt),yes yes)
TEST_SH += \
		$(DIR)/test_NUHMSSMaltEFTHiggs.sh
endif

ifeq ($(WITH_HSSUSY) $(WITH_MSSMEFTHiggs) $(WITH_MSSMMuBMu),yes yes yes)
TEST_SH += \
		$(DIR)/test_MSSMEFTHiggs.sh
endif

ifeq ($(WITH_HSSUSY) $(WITH_MSSMEFTHiggs) $(WITH_NUHMSSMNoFVHimalaya),yes yes yes)
TEST_SH += \
		$(DIR)/test_Mh_uncertainties.sh
endif

ifeq ($(WITH_MSSMEFTHiggs),yes)
TEST_SRC += \
		$(DIR)/test_MSSMEFTHiggs_lambda_threshold_correction.cpp
TEST_SH += \
		$(DIR)/test_MSSMEFTHiggs_librarylink.sh \
		$(DIR)/test_MSSMEFTHiggs_profile.sh
endif

ifeq ($(WITH_MSSMEFTHiggs) $(WITH_MSSMNoFVEFTHiggs),yes yes)
TEST_SH += \
		$(DIR)/test_MSSMNoFVEFTHiggs.sh
endif

ifeq ($(WITH_NMSSMEFTHiggs) $(WITH_lowNMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_NMSSMEFTHiggs.sh
endif

ifeq ($(WITH_SMHighPrecision) $(WITH_SMEFTHiggs),yes yes)
TEST_SH += \
		$(DIR)/test_SMEFTHiggs.sh
endif

ifeq ($(WITH_SM) $(WITH_SMEFTHiggs),yes yes)
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


ifeq ($(WITH_HGTHDMIIMSSMBC) $(WITH_MSSM),yes yes)
TEST_META += \
		$(DIR)/test_HGTHDM_threshold_corrections_scale_invariance.m
endif

ifeq ($(WITH_THDMIIMSSMBC) $(WITH_MSSM),yes yes)
TEST_META += \
		$(DIR)/test_THDM_threshold_corrections_scale_invariance.m
endif

ifeq ($(WITH_THDMIIMSSMBC) $(WITH_HGTHDMIIMSSMBC),yes yes)
TEST_META += \
		$(DIR)/test_HGTHDM_THDM_threshold_corrections_scale_invariance.m
endif

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

$(DIR)/test_cast_model.x: $(DIR)/test_cast_model.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_ckm.x: $(DIR)/test_ckm.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_lowe.x: $(DIR)/test_lowe.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_betafunction.x: $(DIR)/test_betafunction.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_effective_couplings.x: $(DIR)/test_effective_couplings.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_ewsb_solver.x: $(DIR)/test_ewsb_solver.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_fixed_point_iterator.x: $(DIR)/test_fixed_point_iterator.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_goldstones.x: $(DIR)/test_goldstones.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_gsl_vector.x: $(DIR)/test_gsl_vector.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_linalg2.x: $(DIR)/test_linalg2.o
		$(CXX) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(DIR)/test_MSSM_2L_limits.x: $(DIR)/test_MSSM_2L_limits.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_minimizer.x: $(DIR)/test_minimizer.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_mssm_twoloop_mb.x: $(DIR)/test_mssm_twoloop_mb.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_mssm_twoloop_mt.x: $(DIR)/test_mssm_twoloop_mt.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_numerics.x: $(DIR)/test_numerics.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_problems.x: $(DIR)/test_problems.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_pv.x: $(DIR)/test_pv.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_QedQcd.x: $(DIR)/test_QedQcd.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_rk.x: $(DIR)/test_rk.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_root_finder.x: $(DIR)/test_root_finder.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_sminput.x: $(DIR)/test_sminput.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_slha_io.x: $(DIR)/test_slha_io.o $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_thread_pool.x: $(DIR)/test_thread_pool.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_threshold_corrections.x: $(DIR)/test_threshold_corrections.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_threshold_loop_functions.x: $(DIR)/test_threshold_loop_functions.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_spectrum_generator_settings.x: $(DIR)/test_spectrum_generator_settings.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_wrappers.x: $(DIR)/test_wrappers.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_legacy_diagonalization.x: $(DIR)/test_legacy_diagonalization.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_sum.x: $(DIR)/test_sum.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(LIBTEST)

$(DIR)/test_two_scale_mssm_solver.x: $(DIR)/test_two_scale_mssm_solver.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS)

$(DIR)/test_two_scale_mssm_initial_guesser.x: $(DIR)/test_two_scale_mssm_initial_guesser.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS)

$(DIR)/test_two_scale_running_precision.x: $(DIR)/test_two_scale_running_precision.o $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_CMSSM_NMSSM_linking.x: $(DIR)/test_CMSSM_NMSSM_linking.o $(LIBCMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ifeq ($(ENABLE_LOOPTOOLS),yes)
$(DIR)/test_pv_fflite.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(LIBFFLITE)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_looptools.x: $(DIR)/test_pv_crosschecks.cpp $(LIBFLEXI)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(LOOPTOOLSLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_softsusy.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(filter-out %pv.o,$(LIBFLEXI_OBJ))
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)
endif

$(DIR)/test_CMSSMNoFV_benchmark.x: CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)
$(DIR)/test_CMSSMNoFV_benchmark.x: $(DIR)/test_CMSSMNoFV_benchmark.cpp $(RUN_CMSSM_EXE) $(RUN_SOFTPOINT_EXE) $(LIBTEST)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$<) $(LIBTEST) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_compare_ewsb_solvers.x: $(LIBCMSSMGSLHybrid) $(LIBCMSSMGSLHybridS) $(LIBCMSSMGSLBroyden) $(LIBCMSSMGSLNewton) $(LIBCMSSMFPIRelative) $(LIBCMSSMFPIAbsolute) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_loopfunctions.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_sfermions.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_database.x: $(DIR)/test_CMSSM_database.o $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(THREADLIBS)

$(DIR)/test_CMSSM_gluino.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_MRSSM2_gmm2.x: $(LIBMRSSM2) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_model.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_info.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_two_loop_spectrum.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_beta_function_benchmark.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)

$(DIR)/test_CMSSM_initial_guesser.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_higgs_iteration.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_low_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha_output.x: $(DIR)/test_CMSSM_slha_output.o $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(EXAMPLES_EXE) $(DIR)/test_CMSSM_slha_output.in.spc $(RUN_SOFTPOINT_EXE)
		$(CXX) -o $@ $(call abspathx,$< $(LIBCMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_CMSSM_slha_input.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_slha.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCKM_high_scale_constraint.x \
$(DIR)/test_CMSSMCKM_low_scale_constraint.x \
$(DIR)/test_CMSSMCKM_tree_level_spectrum.x: \
	CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

$(DIR)/test_CMSSMCKM_high_scale_constraint.x \
$(DIR)/test_CMSSMCKM_low_scale_constraint.x \
$(DIR)/test_CMSSMCKM_tree_level_spectrum.x: \
	$(LIBSoftsusyFlavourMSSM) $(LIBSoftsusyMSSM) $(LIBCMSSMCKM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCKM_spectrum.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_CMSSM_effective_couplings.x: $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_weinberg_angle.x: $(LIBSoftsusyMSSM) $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSM_weinberg_angle_meta.x: $(LIBCMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMMassWInput_spectrum.x: $(LIBSoftsusyMSSM) $(LIBCMSSMMassWInput) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMLowPrecision.x: $(LIBSoftsusyMSSM) $(LIBCMSSMLowPrecision) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPV_ewsb.x: $(LIBCMSSMCPV) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMCPV) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_MSSMEFTHiggs_lambda_threshold_correction.x: $(LIBMSSMEFTHiggs) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSMCPV_ewsb.x: $(LIBNMSSMCPV) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSMCPV_tree_level_spectrum.x: $(LIBNMSSM) $(LIBNMSSMCPV) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_beta_functions.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBTEST)

$(DIR)/test_NMSSM_ewsb.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_high_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_initial_guesser.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_low_scale_constraint.x: $(LIBSoftsusyNMSSM) $(LIBSoftsusyMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_one_loop_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_self_energies.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_susy_scale_constraint.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_tree_level_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NMSSM_benchmark.x: CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)
$(DIR)/test_NMSSM_benchmark.x: $(DIR)/test_NMSSM_benchmark.cpp $(RUN_NMSSM_EXE) $(RUN_SOFTPOINT_EXE) $(LIBTEST)
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$<) $(LIBTEST) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_NMSSM_slha_output.x: $(DIR)/test_NMSSM_slha_output.o $(LIBNMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(EXAMPLES_EXE) $(DIR)/test_NMSSM_slha_output.in.spc $(RUN_SOFTPOINT_EXE)
		$(CXX) -o $@ $(call abspathx,$< $(LIBNMSSM) $(LIBSoftsusyMSSM) $(LIBFLEXI) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS)

$(DIR)/test_SMSSM_beta_functions.x: $(LIBSMSSM) $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_ewsb.x: $(LIBSMSSM) $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_one_loop_spectrum.x: $(LIBSMSSM) $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSSM_tree_level_spectrum.x: $(LIBSMSSM) $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NUTNMSSM_spectrum.x: $(LIBSoftsusyMSSM) $(LIBSoftsusyNMSSM) $(LIBNUTNMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMNoFV_beta_functions.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_CMSSMNoFV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_CMSSMNoFV_two_loop_spectrum.x: $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_CMSSMNoFV_low_scale_constraint.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_gm2calc.x: $(LIBMSSMNoFVSLHA2) $(LIBGM2Calc) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_MSSMNoFV_onshell.x: $(LIBGM2Calc) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_SM_beta_functions.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_effective_couplings.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_gmm2.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_higgs_loop_corrections.x: $(DIR)/test_SM_higgs_loop_corrections.o $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$< $(LIBSM) $(LIBFLEXI) $(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_SM_low_scale_constraint.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_tree_level_spectrum.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_one_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_three_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_two_loop_spectrum.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_weinberg_angle.x: $(LIBSoftsusyMSSM) $(LIBSM) $(LIBFLEXI) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SM_weinberg_angle_meta.x: $(LIBSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMNoFV_weinberg_angle_meta.x: $(LIBCMSSM) $(LIBCMSSMNoFV) $(LIBFLEXI) $(LIBTEST)

$(DIR)/test_SMHighPrecision_two_loop_spectrum.x: $(LIBSMHighPrecision) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSU3_low_scale_constraint.x: $(LIBSMSU3) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_NSM_low_scale_constraint.x: $(LIBNSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_VCMSSM_ewsb.x: $(LIBVCMSSM) $(LIBCMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMSemiAnalytic_ewsb.x: $(LIBCMSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMSemiAnalytic_consistent_solutions.x: $(LIBCMSSMSemiAnalytic) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMSemiAnalytic_ewsb_solution.x: $(LIBCMSSMSemiAnalytic) $(LIBCMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBCMSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPVSemiAnalytic_ewsb.x: $(LIBCMSSMCPVSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CMSSMCPVSemiAnalytic_semi_analytic_solutions.x: $(LIBCMSSMCPVSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CNMSSM_ewsb.x: $(LIBCNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CNMSSM_semi_analytic_solutions.x: $(LIBCNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CNMSSM_consistent_solutions.x: $(LIBCNMSSM) $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CE6SSM_ewsb.x: $(LIBCE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CE6SSM_semi_analytic_solutions.x: $(LIBCE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_CE6SSM_consistent_solutions.x: $(LIBCE6SSM) $(LIBE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_lowNMSSM_spectrum.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_lowNUHMSSMSemiAnalytic_ewsb.x: $(LIBlowNUHMSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_lowNUHMSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBlowNUHMSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_lowNUHMSSMSemiAnalytic_consistent_solutions.x: $(LIBlowNUHMSSMSemiAnalytic) $(LIBlowNUHMSSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_munuSSMSemiAnalytic_ewsb.x: $(LIBmunuSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_munuSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBmunuSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSemiAnalytic_ewsb.x: $(LIBSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSemiAnalytic_semi_analytic_solutions.x: $(LIBSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SMSemiAnalytic_consistent_solutions.x: $(LIBSMSemiAnalytic) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SSMSemiAnalytic_ewsb.x: $(LIBSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SSMSemiAnalytic_semi_analytic_solutions.x: $(LIBSSMSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_SSMSemiAnalytic_consistent_solutions.x: $(LIBSSMSemiAnalytic) $(LIBSSM) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_ewsb.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_consistent_solutions.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic) $(LIBTHDMII) $(LIBFLEXI) $(LIBLEGACY) $(LIBTEST) $(filter-out -%,$(LOOPFUNCLIBS))

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
