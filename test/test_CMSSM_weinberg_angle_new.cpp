#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle_new

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include <Eigen/Core>

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

void ensure_tree_level_ewsb(CMSSM<Two_scale>& m)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu = m.get_vu();
   const double vd = m.get_vd();
   const double gY = m.get_g1() * sqrt(0.6);
   const double g2 = m.get_g2();
   const double Mu = m.get_Mu();
   const double BMu = m.get_BMu();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
}

// template version of this function already included
// in test_CMSSMNoFV.hpp and test_CMSSM_two_loop_spectrum.cpp
void setup_CMSSM_const(CMSSM<Two_scale>& m, const CMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);
   m.set_thresholds(3);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
}

double CMSSM_low_scale_constraint<Two_scale>::calculate_theta_w(double alpha_em_drbar)
{
   double theta_w = 0.;

   using namespace weinberg_angle;

   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSv = MODELPARAMETER(MSv);

   const double scale         = MODEL->get_scale();
   const double mw_pole       = oneset.displayPoleMW();
   const double mz_pole       = oneset.displayPoleMZ();
   const double mt_pole       = oneset.displayPoleMt();
   const double mt_drbar      = MODEL->get_MFu(2);
   const double mb_drbar      = MODEL->get_MFd(2);
   const double mh_drbar      = MODEL->get_Mhh(0);
   const double gY            = MODEL->get_g1() * 0.7745966692414834;
   const double g2            = MODEL->get_g2();
   const double g3            = MODEL->get_g3();
   const double ymu           = Re(MODEL->get_Ye(1,1));
   const double hmix_12       = MODEL->get_ZH(0,1);
   const double tanBeta       = MODEL->get_vu() / MODEL->get_vd();
   const double mselL         = AbsSqr(ZE(0,0))*MSe(0) + AbsSqr(ZE(1,0))*MSe(1)
      + AbsSqr(ZE(2,0))*MSe(2) + AbsSqr(ZE(3,0))*MSe(3) + AbsSqr(ZE(4,0))*MSe(4)
      + AbsSqr(ZE(5,0))*MSe(5);
   const double msmuL         = AbsSqr(ZE(0,1))*MSe(0) + AbsSqr(ZE(1,1))*MSe(1)
      + AbsSqr(ZE(2,1))*MSe(2) + AbsSqr(ZE(3,1))*MSe(3) + AbsSqr(ZE(4,1))*MSe(4)
      + AbsSqr(ZE(5,1))*MSe(5);
   const double msnue         = AbsSqr(ZV(0,0))*MSv(0) + AbsSqr(ZV(1,0))*MSv(1)
      + AbsSqr(ZV(2,0))*MSv(2);
   const double msnumu        = AbsSqr(ZV(0,1))*MSv(0) + AbsSqr(ZV(1,1))*MSv(1)
      + AbsSqr(ZV(2,1))*MSv(2);
   const double pizztMZ       = Re(MODEL->self_energy_VZ(mz_pole));
   const double piwwt0        = Re(MODEL->self_energy_VWm(0.));
   self_energy_w_at_mw        = Re(MODEL->self_energy_VWm(mw_pole));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY       = gY;
   se_data.g2       = g2;

   double pizztMZ_corrected = pizztMZ;
   double piwwtMW_corrected = self_energy_w_at_mw;
   double piwwt0_corrected  = piwwt0;

   if (model->get_thresholds() > 1) {
      pizztMZ_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole,
            se_data);
      piwwtMW_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(
            self_energy_w_at_mw, mw_pole, se_data);
      piwwt0_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0.,
            se_data);
   }

   Weinberg_angle::Data data;
   data.scale               = scale;
   data.alpha_em_drbar      = ALPHA_EM_DRBAR;
   data.fermi_contant       = oneset.displayFermiConstant();
   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
   data.mw_pole             = mw_pole;
   data.mz_pole             = mz_pole;
   data.mt_pole             = mt_pole;
   data.mh_drbar            = mh_drbar;
   data.hmix_12             = hmix_12;
   data.msel_drbar          = mselL;
   data.msmul_drbar         = msmuL;
   data.msve_drbar          = msnue;
   data.msvm_drbar          = msnumu;
   data.mn_drbar            = MODEL->get_MChi();
   data.mc_drbar            = MODEL->get_MCha();
   data.zn                  = MODEL->get_ZN();
   data.um                  = MODEL->get_UM();
   data.up                  = MODEL->get_UP();
   data.gY                  = gY;
   data.g2                  = g2;
   data.g3                  = g3;
   data.tan_beta            = tanBeta;
   data.ymu                 = ymu;

   Weinberg_angle weinberg;
   weinberg.enable_susy_contributions();
   weinberg.set_number_of_loops(MODEL->get_thresholds());
   weinberg.set_data(data);

   const int error = weinberg.calculate();

   THETAW = ArcSin(weinberg.get_sin_theta());

   if (error)
      MODEL->get_problems().flag_no_rho_convergence();
   else
      MODEL->get_problems().unflag_no_rho_convergence();


   return theta_w;
}

BOOST_AUTO_TEST_CASE( test_1 )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const(model, input);
}
