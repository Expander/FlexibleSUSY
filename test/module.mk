include test/SOFTSUSY/module.mk

DIR      := test
MODNAME  := test
WITH_$(MODNAME) := yes
MODtest_MOD := SM SM_thresholds SplitMSSM MSSM_higgs MSSM_thresholds NMSSM_higgs
MODtest_DEP := $(patsubst %,model_specific/%,$(MODtest_MOD))
MODtest_INC := $(patsubst %,-Imodel_specific/%,$(MODtest_MOD))
MODtest_LIB := $(foreach M,$(MODtest_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

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
		$(DIR)/test_error.cpp \
		$(DIR)/test_fixed_point_iterator.cpp \
		$(DIR)/test_goldstones.cpp \
		$(DIR)/test_gsl_vector.cpp \
		$(DIR)/test_minimizer.cpp \
		$(DIR)/test_namespace_collisions.cpp \
		$(DIR)/test_mssm_twoloop_mb.cpp \
		$(DIR)/test_mssm_twoloop_mt.cpp \
		$(DIR)/test_MSSM_2L_limits.cpp \
		$(DIR)/test_numerics.cpp \
		$(DIR)/test_pmns.cpp \
		$(DIR)/test_problems.cpp \
		$(DIR)/test_pv.cpp \
		$(DIR)/test_raii.cpp \
		$(DIR)/test_root_finder.cpp \
		$(DIR)/test_scan.cpp \
		$(DIR)/test_sm_fourloop_as.cpp \
		$(DIR)/test_sminput.cpp \
		$(DIR)/test_slha_io.cpp \
		$(DIR)/test_sum.cpp \
		$(DIR)/test_threshold_corrections.cpp \
		$(DIR)/test_threshold_loop_functions.cpp \
		$(DIR)/test_spectrum_generator_settings.cpp \
		$(DIR)/test_which.cpp \
		$(DIR)/test_wrappers.cpp

TEST_SH := \
		$(DIR)/test_depgen.sh \
		$(DIR)/test_run_examples.sh \
		$(DIR)/test_run_all_spectrum_generators.sh \
		$(DIR)/test_space_dir.sh \
		$(DIR)/test_wolframscript.sh

TEST_META := \
		$(DIR)/test_BetaFunction.m \
		$(DIR)/test_CConversion.m \
		$(DIR)/test_Constraint.m \
		$(DIR)/test_EWSB.m \
		$(DIR)/test_FunctionModifiers.m \
		$(DIR)/test_LoopFunctions.m \
		$(DIR)/test_MSSM_2L_analytic.m \
		$(DIR)/test_MSSM_2L_mt.m \
		$(DIR)/test_MSSM_2L_yb_softsusy.m \
		$(DIR)/test_MSSM_2L_yt.m \
		$(DIR)/test_MSSM_2L_yt_loopfunction.m \
		$(DIR)/test_MSSM_2L_yt_softsusy.m \
		$(DIR)/test_MRSSM_TreeMasses.m \
		$(DIR)/test_Parameters.m \
		$(DIR)/test_ReadSLHA.m \
		$(DIR)/test_RGIntegrator.m \
		$(DIR)/test_SelfEnergies.m \
		$(DIR)/test_SemiAnalytic.m \
		$(DIR)/test_SM_4loop_as.m \
		$(DIR)/test_SM_higgs_loop_corrections.m \
		$(DIR)/test_SM_higgs_loop_corrections_atab.m \
		$(DIR)/test_TextFormatting.m \
		$(DIR)/test_THDM_threshold_corrections.m \
		$(DIR)/test_THDM_threshold_corrections_gauge.m \
		$(DIR)/test_ThreeLoopQCD.m \
		$(DIR)/test_ThresholdCorrections.m \
		$(DIR)/test_TreeMasses.m \
		$(DIR)/test_TwoLoopNonQCD.m \
		$(DIR)/test_Utils.m \
		$(DIR)/test_Vertices.m \
		$(DIR)/test_Vertices_SortCp.m \
		$(DIR)/test_Vertices_colorsum.m


ifeq ($(ENABLE_THREADS),yes)
TEST_SRC += \
		$(DIR)/test_thread_pool.cpp
endif

ifeq ($(ENABLE_TSIL),yes)
TEST_SRC += \
		$(DIR)/test_sm_twoloop_mt.cpp
endif

ifneq ($(findstring two_scale,$(SOLVERS)),)
TEST_SRC += \
		$(DIR)/test_two_scale_running_precision.cpp \
		$(DIR)/test_two_scale_solver.cpp

ifeq ($(WITH_SOFTSUSY),yes)
TEST_SRC += \
		$(DIR)/test_betafunction.cpp \
		$(DIR)/test_legacy_diagonalization.cpp \
		$(DIR)/test_lowe.cpp \
		$(DIR)/test_QedQcd.cpp \
		$(DIR)/test_rk.cpp \
		$(DIR)/test_two_scale_mssm_solver.cpp \
		$(DIR)/test_two_scale_mssm_initial_guesser.cpp
endif

ifeq ($(WITH_SM) $(WITH_SOFTSUSY),yes yes)
TEST_SRC += \
		$(DIR)/test_SM_weinberg_angle.cpp \
		$(DIR)/test_SM_weinberg_angle_meta.cpp
endif

ifeq ($(WITH_SOFTSUSY) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_loopfunctions.cpp \
		$(DIR)/test_sfermions.cpp \
		$(DIR)/test_CMSSM_beta_function_benchmark.cpp \
		$(DIR)/test_CMSSM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSM_higgs_iteration.cpp \
		$(DIR)/test_CMSSM_initial_guesser.cpp \
		$(DIR)/test_CMSSM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSM_model.cpp \
		$(DIR)/test_CMSSM_slha_output.cpp \
		$(DIR)/test_CMSSM_spectrum.cpp \
		$(DIR)/test_CMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_CMSSM_weinberg_angle.cpp \
		$(DIR)/test_CMSSM_weinberg_angle_meta.cpp
TEST_SH += \
		$(DIR)/test_CMSSM_gluino.sh
endif

ifeq ($(WITH_SOFTSUSY) $(WITH_CMSSMMassWInput),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMMassWInput_spectrum.cpp
endif

ifeq ($(WITH_CMSSMLowPrecision),yes)
TEST_SRC += \
		$(DIR)/test_CMSSMLowPrecision.cpp
endif

ifeq ($(WITH_SOFTSUSY) $(WITH_CMSSMCKM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCKM_high_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_low_scale_constraint.cpp \
		$(DIR)/test_CMSSMCKM_tree_level_spectrum.cpp
TEST_SH += \
		$(DIR)/test_CMSSMCKM_spectrum.sh
endif

ifeq ($(WITH_SOFTSUSY) $(WITH_NMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_NMSSM_benchmark.cpp \
		$(DIR)/test_NMSSM_beta_functions.cpp \
		$(DIR)/test_NMSSM_ewsb.cpp \
		$(DIR)/test_NMSSM_high_scale_constraint.cpp \
		$(DIR)/test_NMSSM_initial_guesser.cpp \
		$(DIR)/test_NMSSM_low_scale_constraint.cpp \
		$(DIR)/test_NMSSM_one_loop_spectrum.cpp \
		$(DIR)/test_NMSSM_self_energies.cpp \
		$(DIR)/test_NMSSM_slha_output.cpp \
		$(DIR)/test_NMSSM_spectrum.cpp \
		$(DIR)/test_NMSSM_susy_scale_constraint.cpp \
		$(DIR)/test_NMSSM_tree_level_spectrum.cpp
endif

ifeq ($(WITH_SOFTSUSY) $(WITH_SMSSM),yes yes)
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

ifeq ($(WITH_CMSSMCKM),yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCKM_b_to_s_gamma_internal_spectrum.cpp
endif

ifeq ($(WITH_MRSSM2),yes)
TEST_SRC += \
		$(DIR)/test_MRSSM2_gmm2.cpp \
		$(DIR)/test_MRSSM2_mw_calculation.cpp \
		$(DIR)/test_MRSSM2_l_to_lgamma.cpp
endif

ifeq ($(WITH_MRSSM2CKM),yes)
TEST_SRC += \
		$(DIR)/test_MRSSM2CKM_b_to_s_gamma.cpp
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

ifeq ($(WITH_munuSSM), yes)
TEST_META += \
		$(DIR)/test_munuSSM_TreeMasses.m

TEST_SRC += \
		$(DIR)/test_munuSSM_gmm2.cpp
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

ifeq ($(WITH_SOFTSUSY) $(WITH_CMSSMNoFV),yes yes)
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

ifeq ($(WITH_SOFTSUSY) $(WITH_CMSSM) $(WITH_CMSSMNoFV),yes yes yes)
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

$(DIR)/test_pv_crosschecks.sh.xml: $(TEST_PV_EXE)

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
		$(DIR)/test_CMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/test_CMSSM_slha.cpp \
		$(DIR)/test_CMSSM_slha_input.cpp \
		$(DIR)/test_CMSSM_two_loop_spectrum.cpp \
		$(DIR)/test_CMSSM_info.cpp \
		$(DIR)/test_CMSSM_mw_calculation.cpp
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
		$(DIR)/test_SM_mass_eigenstates_interface.cpp \
		$(DIR)/test_SM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/test_SM_one_loop_spectrum.cpp \
		$(DIR)/test_SM_higgs_loop_corrections.cpp \
		$(DIR)/test_SM_tree_level_spectrum.cpp \
		$(DIR)/test_SM_three_loop_spectrum.cpp \
		$(DIR)/test_SM_two_loop_spectrum.cpp \
		$(DIR)/test_SM_mw_calculation.cpp
#		$(DIR)/test_SM_cxxdiagrams.cpp
endif

ifeq ($(ENABLE_FEYNARTS) $(ENABLE_FORMCALC),yes yes)
ifeq ($(WITH_SM),yes)
# TEST_SRC += \
#		$(DIR)/test_SM_npointfunctions.cpp \
#		$(DIR)/test_SM_matching_selfenergy_Fd.cpp
endif
ifeq ($(WITH_MSSM),yes)
#TEST_SRC += \
#		$(DIR)/test_MSSM_npointfunctions.cpp \
#		$(DIR)/test_MSSM_matching_selfenergy_Fd.cpp
endif
endif

ifeq ($(WITH_SMHighPrecision),yes)
TEST_SRC += \
		$(DIR)/test_SMHighPrecision_two_loop_spectrum.cpp
endif

ifeq ($(WITH_SMSU3),yes)
TEST_SRC += \
		$(DIR)/test_SMSU3_low_scale_constraint.cpp
TEST_META += \
		$(DIR)/test_SMSU3_TreeMasses.m
endif

ifeq ($(WITH_NSM),yes)
TEST_SRC += \
		$(DIR)/test_NSM_low_scale_constraint.cpp
endif

ifeq ($(WITH_lowMSSM) $(WITH_CMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_lowMSSM.sh
endif

ifeq ($(WITH_lowNMSSM) $(WITH_SOFTSUSY),yes yes)
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
ifeq ($(WITH_CMSSMCPV) $(ENABLE_LIBRARYLINK),yes yes)
TEST_META += \
		$(DIR)/test_CMSSMCPV_librarylink.m
endif

ifeq ($(WITH_CMSSMCPV) $(WITH_CMSSM),yes yes)
TEST_SRC += \
		$(DIR)/test_CMSSMCPV_tree_level_spectrum.cpp
TEST_SH += \
		$(DIR)/test_CMSSMCPV_spectrum.sh
endif

ifeq ($(WITH_MSSMNoFV),yes)
TEST_META += \
		$(DIR)/test_MSSMNoFV_TreeMasses.m
endif

ifeq ($(WITH_MSSMCPV),yes)
TEST_META += \
		$(DIR)/test_MSSMCPV_TreeMasses.m
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

ifeq ($(WITH_NUTNMSSM) $(WITH_SOFTSUSY),yes yes)
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
TEST_META += \
		$(DIR)/test_HSSUSY_uncertainty.m
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
TEST_META += \
		$(DIR)/test_MSSMEFTHiggs_uncertainty.m
TEST_SRC += \
		$(DIR)/test_MSSMEFTHiggs_lambda_threshold_correction.cpp
TEST_SH += \
		$(DIR)/test_MSSMEFTHiggs_profile.sh
endif
ifeq ($(WITH_MSSMEFTHiggs) $(ENABLE_LIBRARYLINK),yes yes)
TEST_SH += \
		$(DIR)/test_MSSMEFTHiggs_librarylink.sh
endif

ifeq ($(WITH_MSSMEFTHiggs) $(WITH_MSSMNoFVEFTHiggs),yes yes)
TEST_SH += \
		$(DIR)/test_MSSMNoFVEFTHiggs.sh
endif

ifeq ($(WITH_NMSSMEFTHiggs) $(WITH_lowNMSSM),yes yes)
TEST_SH += \
		$(DIR)/test_NMSSMEFTHiggs.sh
endif

ifeq ($(WITH_NUHMSSMNoFVHimalaya),yes)
TEST_META += \
		$(DIR)/test_NUHMSSMNoFVHimalaya_uncertainty.m
endif

ifeq ($(WITH_SMHighPrecision) $(WITH_SMEFTHiggs),yes yes)
TEST_SH += \
		$(DIR)/test_SMEFTHiggs.sh
endif

ifeq ($(WITH_SplitMSSMEFTHiggs),yes)
TEST_SH += \
		$(DIR)/test_SplitMSSMEFTHiggs.sh
endif

ifeq ($(WITH_SM) $(WITH_SMEFTHiggs) $(ENABLE_LIBRARYLINK),yes yes yes)
TEST_META += \
		$(DIR)/test_multiple_librarylinks.m
endif

ifeq ($(WITH_SM) $(ENABLE_LIBRARYLINK),yes yes)
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
		$(DIR)/test_SM_3loop_beta.m \
		$(DIR)/test_SM_TreeMasses.m
endif

ifeq ($(WITH_SMnoGUT),yes)
TEST_META += \
		$(DIR)/test_SMnoGUT_3loop_beta.m
endif

ifeq ($(WITH_CMSSM),yes)
TEST_META += \
		$(DIR)/test_CMSSM_3loop_beta.m
endif
ifeq ($(WITH_CMSSM) $(ENABLE_LIBRARYLINK),yes yes)
TEST_META += \
		$(DIR)/test_CMSSM_librarylink.m \
		$(DIR)/test_CMSSM_librarylink_parallel.m
TEST_SH += \
		$(DIR)/test_CMSSM_librarylink_slha.sh
endif

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_OBJ +=	$(foreach i, 1 2 3 4 5 6 7 8, $(DIR)/test_linalg2_part$(i).o)

$(DIR)/test_linalg2_part%.d: $(DIR)/test_linalg2.cpp | $(DEPGEN)
	$(Q)$(DEPGEN) $(CPPFLAGS) -DTEST_LINALG2_PART$* -MM -MI -o '$@' -MT '$*.o' $^

$(DIR)/test_linalg2_part%.o: $(DIR)/test_linalg2.cpp
	@$(MSG)
	$(Q)$(CXX) $(CPPFLAGS) -DTEST_LINALG2_PART$* $(CXXFLAGS) -c $< -o $@

$(foreach i, 2 3 4 5 6 7 8, $(eval \
$(DIR)/test_linalg2_part$(i).o: | $(DIR)/test_linalg2_part$(shell echo $$(($(i)-1))).o))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_XML      := $(DIR)/test.xml
TEST_EXE_XML  := $(TEST_EXE:.x=.x.xml)
TEST_SH_XML   := $(TEST_SH:.sh=.sh.xml)
TEST_META_XML := $(TEST_META:.m=.m.xml)
TEST_ALL_XML  := $(TEST_EXE_XML) $(TEST_SH_XML)
ifeq ($(ENABLE_META),yes)
TEST_ALL_XML  += $(TEST_META_XML)
endif

TEST_ALL_LOG  := $(TEST_ALL_XML:.xml=.log)

ifeq ($(ENABLE_LOOPTOOLS),yes)
TEST_EXE += $(TEST_PV_EXE)

$(DIR)/test_pv_fflite.x    : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_FFLITE
$(DIR)/test_pv_looptools.x : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_LOOPTOOLS
$(DIR)/test_pv_softsusy.x  : CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS) -DTEST_PV_SOFTSUSY
endif

$(DIR)/test_threshold_loop_functions.x: CPPFLAGS += -DTEST_DATA_DIR="\"test/data/threshold_loop_functions\""

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-log \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj \
		execute-tests execute-meta-tests execute-compiled-tests \
		execute-shell-tests

all-$(MODNAME): $(LIBTEST) $(TEST_EXE) $(TEST_XML)
		@printf "%s\n" "All tests passed."

clean-$(MODNAME)-dep: clean-SOFTSUSY-dep
		$(Q)-rm -f $(TEST_DEP)
		$(Q)-rm -f $(LIBTEST_DEP)

clean-$(MODNAME)-lib: clean-SOFTSUSY-lib
		$(Q)-rm -f $(LIBTEST)

clean-$(MODNAME)-obj: clean-SOFTSUSY-obj
		$(Q)-rm -f $(TEST_OBJ)
		$(Q)-rm -f $(LIBTEST_OBJ)

clean-$(MODNAME)-log:
		$(Q)-rm -f $(TEST_XML)
		$(Q)-rm -f $(TEST_ALL_XML)
		$(Q)-rm -f $(TEST_ALL_LOG)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
                  clean-$(MODNAME)-lib clean-$(MODNAME)-log
		$(Q)-rm -f $(TEST_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

execute-tests:  $(TEST_XML)
		@printf "%s\n" "All tests passed."

ifeq ($(ENABLE_META),yes)
execute-meta-tests: $(TEST_META_XML)
		@printf "%s\n" "All meta tests passed."
else
execute-meta-tests:
		@printf "%s\n" "All meta tests passed."
endif

execute-compiled-tests: $(TEST_EXE_XML)
		@printf "%s\n" "All compiled tests passed."

execute-shell-tests: $(TEST_SH_XML)
		@printf "%s\n" "All shell script tests passed."

# creates .xml file with test result
PTR = write_test_result_file() { \
	echo "\
<test>\n\
\t<name>$$2</name>\n\
\t<date>$$(date)</date>\n\
\t<commit>$$(git describe --tags 2> /dev/null || echo unknown)</commit>\n\
\t<status>$$1</status>\n\
\t<logfile>$$(basename $$4)</logfile>\n\
</test>" > "$$3" ; \
	if [ $$1 = 0 ]; then \
		printf "%-66s %4s\n" "$$2" "OK"; \
	else \
		printf "%-66s %4s\n" "$$2" "FAILED"; \
	fi; \
	return $$1; \
}

$(DIR)/%.x.xml: $(DIR)/%.x
		@rm -f $@ $(@:.xml=.log)
		@$(PTR); \
		BOOST_TEST_CATCH_SYSTEM_ERRORS="no" \
		$< --log_level=test_suite >> $(@:.xml=.log) 2>&1; \
		write_test_result_file $$? $< $@ $(@:.xml=.log)

$(DIR)/%.m.xml: $(DIR)/%.m $(META_SRC)
		@rm -f $@ $(@:.xml=.log)
		@$(PTR); \
		printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; \
		Quit[TestSuite\`GetNumberOfFailedTests[]]" | "$(MATH)" >> $(@:.xml=.log) 2>&1; \
		write_test_result_file $$? $< $@ $(@:.xml=.log)

$(DIR)/%.sh.xml: $(DIR)/%.sh
		@rm -f $@ $(@:.xml=.log)
		@$(PTR); \
		MATH_CMD="$(MATH)" $< >> $(@:.xml=.log) 2>&1; \
		write_test_result_file $$? $< $@ $(@:.xml=.log)

$(TEST_XML): $(TEST_ALL_XML)
		@echo "Creating $@"
		@echo "\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<?xml-stylesheet type=\"text/xsl\" href=\"test.xsl\"?>\n\
<tests date=\"$$(date)\">\n\
$$(for f in $^ ; do echo "\t<test filename=\"$$(basename $$f)\"/>"; done)\n\
</tests>" > $@

$(DIR)/test_depgen.sh.xml: $(DEPGEN_EXE)

$(DIR)/test_lowMSSM.sh.xml: $(RUN_CMSSM_EXE) $(RUN_lowMSSM_EXE)

$(DIR)/test_run_all_spectrum_generators.sh.xml: allexec

$(DIR)/test_CMSSM_NMSSM_linking.x: $(LIBCMSSM) $(LIBNMSSM)

ifeq ($(ENABLE_LOOPTOOLS),yes)
$(DIR)/test_pv_fflite.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(LIBFFLITE)
		@$(MSG)
		$(Q)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_looptools.x: $(DIR)/test_pv_crosschecks.cpp $(LIBFLEXI)
		@$(MSG)
		$(Q)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(LOOPFUNCLIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(FLIBS)

$(DIR)/test_pv_softsusy.x: $(DIR)/test_pv_crosschecks.cpp src/pv.cpp $(filter-out %pv.o,$(LIBFLEXI_OBJ))
		@$(MSG)
		$(Q)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(call abspathx,$^) $(SQLITELIBS) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(FLIBS)
endif

$(DIR)/test_CMSSMNoFV_benchmark.x.xml: $(RUN_CMSSM_EXE) $(RUN_SOFTPOINT_EXE)

$(DIR)/test_compare_ewsb_solvers.x: \
	$(LIBCMSSMGSLHybrid) $(LIBCMSSMGSLHybridS) $(LIBCMSSMGSLBroyden) \
	$(LIBCMSSMGSLNewton) $(LIBCMSSMFPIRelative) $(LIBCMSSMFPIAbsolute) \
	$(LIBCMSSMFPITadpole)

$(DIR)/test_loopfunctions.x: $(LIBCMSSM)

$(DIR)/test_sfermions.x: $(LIBCMSSM)

# $(DIR)/test_SM_cxxdiagrams.cpp : $(DIR)/test_SM_cxxdiagrams.meta $(DIR)/test_SM_cxxdiagrams.cpp.in $(META_SRC) $(METACODE_STAMP_SM)
#		@$(MSG)
#		$(Q)printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; Quit[0]" | "$(MATH)"
# $(DIR)/test_SM_cxxdiagrams.x: $(LIBSM)

# $(DIR)/test_SM_npointfunctions.cpp : $(DIR)/test_SM_npointfunctions.meta $(DIR)/test_SM_npointfunctions.cpp.in $(META_SRC) $(METACODE_STAMP_SM)
#		@$(MSG)
#		$(Q)printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; Quit[0]" | "$(MATH)"
# $(DIR)/test_SM_npointfunctions.x: $(LIBSM)

# $(DIR)/test_SM_matching_selfenergy_Fd.cpp : $(DIR)/test_SM_matching_selfenergy_Fd.meta $(DIR)/test_SM_matching_selfenergy_Fd.cpp.in $(META_SRC) $(METACODE_STAMP_SM)
# 		printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; Quit[0]" | "$(MATH)"
# $(DIR)/test_SM_matching_selfenergy_Fd.x: $(LIBSM)

# $(DIR)/test_MSSM_npointfunctions.cpp : \
#		$(DIR)/test_MSSM_npointfunctions.meta \
#		$(DIR)/test_MSSM_npointfunctions.cpp.in \
#		$(META_SRC) $(METACODE_STAMP_MSSM)
#		@$(MSG)
#		$(Q)printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; Quit[0]" | "$(MATH)"
# $(DIR)/test_MSSM_npointfunctions.x: $(LIBMSSM)

# $(DIR)/test_MSSM_matching_selfenergy_Fd.cpp : \
#		$(DIR)/test_MSSM_matching_selfenergy_Fd.meta \
#		$(DIR)/test_MSSM_matching_selfenergy_Fd.cpp.in \
#		$(META_SRC) $(METACODE_STAMP_MSSM)
#		@$(MSG)
#		$(Q)printf "%s" "AppendTo[\$$Path, \"./meta/\"]; Get[\"$<\"]; Quit[0]" | "$(MATH)"
# $(DIR)/test_MSSM_matching_selfenergy_Fd.x: $(LIBMSSM)

$(DIR)/test_CMSSM_database.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_gluino.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_MRSSM2_gmm2.x: $(LIBMRSSM2)

$(DIR)/test_CMSSM_mass_eigenstates_decoupling_scheme.x: $(LIBCMSSM)

$(DIR)/test_MRSSM2_mw_calculation.x: $(LIBMRSSM2)

$(DIR)/test_MRSSM2_l_to_lgamma.x: $(LIBMRSSM2)

$(DIR)/test_MRSSM2CKM_b_to_s_gamma.x: $(LIBMRSSM2CKM)

$(DIR)/test_CMSSMCKM_b_to_s_gamma_internal_spectrum.x: $(LIBCMSSMCKM)

$(DIR)/test_CMSSM_model.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_info.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_two_loop_spectrum.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_mw_calculation.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_beta_function_benchmark.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_initial_guesser.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_higgs_iteration.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_high_scale_constraint.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_low_scale_constraint.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_susy_scale_constraint.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_slha_output.x.xml: $(EXAMPLES_EXE) $(DIR)/test_CMSSM_slha_output.in.spc $(RUN_SOFTPOINT_EXE)
$(DIR)/test_CMSSM_slha_output.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_slha_input.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_slha.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_spectrum.x: $(LIBCMSSM)

$(DIR)/test_CMSSMCKM_high_scale_constraint.x \
$(DIR)/test_CMSSMCKM_low_scale_constraint.x \
$(DIR)/test_CMSSMCKM_tree_level_spectrum.x: \
	 $(LIBCMSSMCKM)

$(DIR)/test_CMSSMCKM_spectrum.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_CMSSM_effective_couplings.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_weinberg_angle.x: $(LIBCMSSM)

$(DIR)/test_CMSSM_weinberg_angle_meta.x: $(LIBCMSSM)

$(DIR)/test_CMSSMMassWInput_spectrum.x: $(LIBCMSSMMassWInput)

$(DIR)/test_CMSSMLowPrecision.x: $(LIBCMSSMLowPrecision)

$(DIR)/test_CMSSMCPV_ewsb.x: $(LIBCMSSMCPV)

$(DIR)/test_CMSSMCPV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMCPV)

$(DIR)/test_MSSMEFTHiggs_lambda_threshold_correction.x: $(LIBMSSMEFTHiggs)

$(DIR)/test_NMSSMCPV_ewsb.x: $(LIBNMSSMCPV)

$(DIR)/test_NMSSMCPV_tree_level_spectrum.x: $(LIBNMSSM) $(LIBNMSSMCPV)

$(DIR)/test_NMSSM_beta_functions.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_ewsb.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_high_scale_constraint.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_initial_guesser.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_low_scale_constraint.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_one_loop_spectrum.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_self_energies.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_spectrum.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_susy_scale_constraint.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_tree_level_spectrum.x: $(LIBNMSSM)

$(DIR)/test_NMSSM_benchmark.x.xml: $(RUN_NMSSM_EXE) $(RUN_SOFTPOINT_EXE)

$(DIR)/test_NMSSM_slha_output.x.xml: $(EXAMPLES_EXE) $(DIR)/test_NMSSM_slha_output.in.spc $(RUN_SOFTPOINT_EXE)
$(DIR)/test_NMSSM_slha_output.x: $(LIBNMSSM)

$(DIR)/test_SMSSM_beta_functions.x: $(LIBSMSSM)

$(DIR)/test_SMSSM_ewsb.x: $(LIBSMSSM)

$(DIR)/test_SMSSM_one_loop_spectrum.x: $(LIBSMSSM)

$(DIR)/test_SMSSM_tree_level_spectrum.x: $(LIBSMSSM)

$(DIR)/test_NUTNMSSM_spectrum.x: $(LIBNUTNMSSM)

$(DIR)/test_CMSSMNoFV_beta_functions.x: $(LIBCMSSM) $(LIBCMSSMNoFV)

$(DIR)/test_CMSSMNoFV_tree_level_spectrum.x: $(LIBCMSSM) $(LIBCMSSMNoFV)

$(DIR)/test_CMSSMNoFV_two_loop_spectrum.x: $(LIBCMSSMNoFV)

$(DIR)/test_CMSSMNoFV_low_scale_constraint.x: $(LIBCMSSM) $(LIBCMSSMNoFV)

$(DIR)/test_gm2calc.x: $(LIBMSSMNoFVSLHA2) $(LIBGM2Calc)

$(DIR)/test_MSSMNoFV_onshell.x: $(LIBGM2Calc)

$(DIR)/test_SM_beta_functions.x: $(LIBSM)

$(DIR)/test_SM_effective_couplings.x: $(LIBSM)

$(DIR)/test_SM_gmm2.x: $(LIBSM)

$(DIR)/test_SM_higgs_loop_corrections.x: $(LIBSM)

$(DIR)/test_SM_low_scale_constraint.x: $(LIBSM)

$(DIR)/test_SM_mass_eigenstates_interface.x: $(LIBSM)

$(DIR)/test_SM_mass_eigenstates_decoupling_scheme.x: $(LIBSM)

$(DIR)/test_SM_tree_level_spectrum.x: $(LIBSM)

$(DIR)/test_SM_one_loop_spectrum.x: $(LIBSM)

$(DIR)/test_SM_three_loop_spectrum.x: $(LIBSM)

$(DIR)/test_SM_two_loop_spectrum.x: $(LIBSM)

$(DIR)/test_SM_mw_calculation.x: $(LIBSM)

$(DIR)/test_SM_weinberg_angle.x: $(LIBSM)

$(DIR)/test_SM_weinberg_angle_meta.x: $(LIBSM)

$(DIR)/test_CMSSMNoFV_weinberg_angle_meta.x: $(LIBCMSSM) $(LIBCMSSMNoFV)

$(DIR)/test_SMHighPrecision_two_loop_spectrum.x: $(LIBSMHighPrecision)

$(DIR)/test_SMSU3_low_scale_constraint.x: $(LIBSMSU3)

$(DIR)/test_NSM_low_scale_constraint.x: $(LIBNSM)

$(DIR)/test_VCMSSM_ewsb.x: $(LIBVCMSSM) $(LIBCMSSM)

$(DIR)/test_CMSSMSemiAnalytic_ewsb.x: $(LIBCMSSMSemiAnalytic)

$(DIR)/test_CMSSMSemiAnalytic_consistent_solutions.x: $(LIBCMSSMSemiAnalytic) $(LIBCMSSM)

$(DIR)/test_CMSSMSemiAnalytic_ewsb_solution.x: $(LIBCMSSMSemiAnalytic) $(LIBCMSSM)

$(DIR)/test_CMSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBCMSSMSemiAnalytic)

$(DIR)/test_CMSSMCPVSemiAnalytic_ewsb.x: $(LIBCMSSMCPVSemiAnalytic)

$(DIR)/test_CMSSMCPVSemiAnalytic_semi_analytic_solutions.x: $(LIBCMSSMCPVSemiAnalytic)

$(DIR)/test_CNMSSM_ewsb.x: $(LIBCNMSSM)

$(DIR)/test_CNMSSM_semi_analytic_solutions.x: $(LIBCNMSSM)

$(DIR)/test_CNMSSM_consistent_solutions.x: $(LIBCNMSSM) $(LIBNMSSM)

$(DIR)/test_CE6SSM_ewsb.x: $(LIBCE6SSM)

$(DIR)/test_CE6SSM_semi_analytic_solutions.x: $(LIBCE6SSM)

$(DIR)/test_CE6SSM_consistent_solutions.x: $(LIBCE6SSM) $(LIBE6SSM)

$(DIR)/test_lowNMSSM_spectrum.sh: $(RUN_SOFTPOINT_EXE)

$(DIR)/test_lowNUHMSSMSemiAnalytic_ewsb.x: $(LIBlowNUHMSSMSemiAnalytic)

$(DIR)/test_lowNUHMSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBlowNUHMSSMSemiAnalytic)

$(DIR)/test_lowNUHMSSMSemiAnalytic_consistent_solutions.x: $(LIBlowNUHMSSMSemiAnalytic) $(LIBlowNUHMSSM)

$(DIR)/test_munuSSM_gmm2.x: $(LIBmunuSSM)

$(DIR)/test_munuSSMSemiAnalytic_ewsb.x: $(LIBmunuSSMSemiAnalytic)

$(DIR)/test_munuSSMSemiAnalytic_semi_analytic_solutions.x: $(LIBmunuSSMSemiAnalytic)

$(DIR)/test_SMSemiAnalytic_ewsb.x: $(LIBSMSemiAnalytic)

$(DIR)/test_SMSemiAnalytic_semi_analytic_solutions.x: $(LIBSMSemiAnalytic)

$(DIR)/test_SMSemiAnalytic_consistent_solutions.x: $(LIBSMSemiAnalytic) $(LIBSM)

$(DIR)/test_SSMSemiAnalytic_ewsb.x: $(LIBSSMSemiAnalytic)

$(DIR)/test_SSMSemiAnalytic_semi_analytic_solutions.x: $(LIBSSMSemiAnalytic)

$(DIR)/test_SSMSemiAnalytic_consistent_solutions.x: $(LIBSSMSemiAnalytic) $(LIBSSM)

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_ewsb.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic)

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic)

$(DIR)/test_THDMIIEWSBAtMZSemiAnalytic_consistent_solutions.x: $(LIBTHDMIIEWSBAtMZSemiAnalytic) $(LIBTHDMII)

# adding libraries to the end of the list of dependencies
$(TEST_EXE): $(LIBSOFTSUSY) $(MODtest_LIB) $(LIBTEST) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))

# general test rule
$(DIR)/test_%.x: $(DIR)/test_%.o
		@$(MSG)
		$(Q)$(CXX) -o $@ $(call abspathx,$^) \
		$(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) \
		$(THREADLIBS) $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS)

# add boost and eigen flags for the test object files and dependencies
$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += -Itest/SOFTSUSY $(MODtest_INC) $(BOOSTFLAGS) $(EIGENFLAGS) $(GSLFLAGS) $(TSILFLAGS)

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIBTEST): $(LIBTEST_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(FLIBS)
else
$(LIBTEST): $(LIBTEST_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIBTEST_DEP) $(TEST_DEP)
ALLLIB += $(LIBTEST)
ALLTST += $(TEST_EXE)
ALLMODDEP += $(MODtest_DEP)
