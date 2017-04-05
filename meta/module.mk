DIR          := meta
MODNAME      := meta

META_MSSM_SRC:= \
		$(DIR)/MSSM/BDSZHiggs.m \
		$(DIR)/MSSM/beta_BMu.m \
		$(DIR)/MSSM/beta_g1.m \
		$(DIR)/MSSM/beta_g2.m \
		$(DIR)/MSSM/beta_g3.m \
		$(DIR)/MSSM/beta_M1.m \
		$(DIR)/MSSM/beta_M2.m \
		$(DIR)/MSSM/beta_M3.m \
		$(DIR)/MSSM/beta_md2.m \
		$(DIR)/MSSM/beta_me2.m \
		$(DIR)/MSSM/beta_mHd2.m \
		$(DIR)/MSSM/beta_mHu2.m \
		$(DIR)/MSSM/beta_ml2.m \
		$(DIR)/MSSM/beta_mq2.m \
		$(DIR)/MSSM/beta_mu2.m \
		$(DIR)/MSSM/beta_Mu.m \
		$(DIR)/MSSM/beta_TYd.m \
		$(DIR)/MSSM/beta_TYe.m \
		$(DIR)/MSSM/beta_TYu.m \
		$(DIR)/MSSM/beta_Yd.m \
		$(DIR)/MSSM/beta_Ye.m \
		$(DIR)/MSSM/beta_Yu.m \
		$(DIR)/MSSM/dmtas2.m \
		$(DIR)/MSSM/extract_MSSM_beta_functions_from_hep-ph-0308231.m \
		$(DIR)/MSSM/extract_MSSM_Mt_over_mt_from_softsusy.m \
		$(DIR)/MSSM/gamma_SdR.m \
		$(DIR)/MSSM/gamma_SeR.m \
		$(DIR)/MSSM/gamma_SHd.m \
		$(DIR)/MSSM/gamma_SHu.m \
		$(DIR)/MSSM/gamma_SlL.m \
		$(DIR)/MSSM/gamma_SqL.m \
		$(DIR)/MSSM/gamma_SuR.m \
		$(DIR)/MSSM/tquark_1loop_qcd.m \
		$(DIR)/MSSM/tquark_1loop_strong.m \
		$(DIR)/MSSM/tquark_2loop_qcd.m \
		$(DIR)/MSSM/tquark_2loop_strong.m \
		$(DIR)/MSSM/tquark_to_cpp.m

META_SM_SRC  := \
		$(DIR)/SM/beta_g1.m \
		$(DIR)/SM/beta_g2.m \
		$(DIR)/SM/beta_g3.m \
		$(DIR)/SM/beta_gb.m \
		$(DIR)/SM/beta_gtau.m \
		$(DIR)/SM/beta_gt.m \
		$(DIR)/SM/beta_lambda.m \
		$(DIR)/SM/beta_m2.m \
		$(DIR)/SM/HSSUSY_corrections.m

META_THDM_SRC:= \
		$(DIR)/THDM/Thresholds_1L_full.m

META_SRC     := \
		$(DIR)/AnomalousDimension.m \
		$(DIR)/BetaFunction.m \
		$(DIR)/CConversion.m \
		$(DIR)/Constraint.m \
		$(DIR)/ConvergenceTester.m \
		$(DIR)/EffectiveCouplings.m \
		$(DIR)/EWSB.m \
		$(DIR)/FlexibleEFTHiggsMatching.m \
		$(DIR)/FlexibleSUSY.m \
		$(DIR)/FlexibleTower.m \
		$(DIR)/Format.m \
		$(DIR)/FSMathLink.m \
		$(DIR)/GMuonMinus2.m \
		$(DIR)/LatticeUtils.m \
		$(DIR)/LoopFunctions.m \
		$(DIR)/LoopMasses.m \
		$(DIR)/Observables.m \
		$(DIR)/Parameters.m \
		$(DIR)/Phases.m \
		$(DIR)/ReadSLHA.m \
		$(DIR)/RGIntegrator.m \
		$(DIR)/SelfEnergies.m \
		$(DIR)/TestSuite.m \
		$(DIR)/TextFormatting.m \
		$(DIR)/ThreeLoopMSSM.m \
		$(DIR)/ThreeLoopQCD.m \
		$(DIR)/ThreeLoopSM.m \
		$(DIR)/ThresholdCorrections.m \
		$(DIR)/Traces.m \
		$(DIR)/TreeMasses.m \
		$(DIR)/TwoLoopSM.m \
		$(DIR)/TwoLoopMSSM.m \
		$(DIR)/TwoLoopQCD.m \
		$(DIR)/Utils.m \
		$(DIR)/Vertices.m \
		$(DIR)/WeinbergAngle.m \
		$(DIR)/WriteOut.m \
		$(DIR)/writeRGE.m \
		$(DIR)/writeNRGE.m \
		$(META_MSSM_SRC) \
		$(META_SM_SRC) \
		$(META_THDM_SRC)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

clean-$(MODNAME):
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
