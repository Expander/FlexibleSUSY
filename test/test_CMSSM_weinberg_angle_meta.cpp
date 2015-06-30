// this test file is almost identical to test_CMSSM_weinberg_angle_pointer.cpp

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle_meta

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "stopwatch.hpp"

#define private public

#include "weinberg_angle.hpp"
#include "CMSSM_weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

void ensure_tree_level_ewsb(CMSSM<Two_scale>& model)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu   = model.get_vu();
   const double vd   = model.get_vd();
   const double gY   = model.get_g1() * Sqrt(0.6);
   const double g2   = model.get_g2();
   const double Mu   = model.get_Mu();
   const double BMu  = model.get_BMu();
   const double mHd2 = BMu*vu/vd - (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   const double mHu2 = BMu*vd/vu + (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);
}

// template version of this function already included
// in test_CMSSMNoFV.hpp and test_CMSSM_two_loop_spectrum.cpp
void setup_CMSSM_const(CMSSM<Two_scale>& model, const CMSSM_input_parameters& input)
{
   const double ALPHASMZ = Electroweak_constants::alpha3;
   const double ALPHAMZ  = Electroweak_constants::aem;
   const double sinthWsq = Electroweak_constants::sinThetaW2;
   const double alpha1   = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2   = ALPHAMZ / sinthWsq;
   const double g1       = Sqrt(4 * Pi * alpha1);
   const double g2       = Sqrt(4 * Pi * alpha2);
   const double g3       = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta  = input.TanBeta;
   const double sinBeta  = sin(atan(tanBeta));
   const double cosBeta  = cos(atan(tanBeta));
   const double M12      = input.m12;
   const double m0       = input.m0;
   const double a0       = input.Azero;
   const double root2    = Sqrt(2.0);
   const double vev      = Electroweak_constants::vev;
   const double vu       = vev * sinBeta;
   const double vd       = vev * cosBeta;
   const double susyMu   = input.SignMu * 120.0;
   const double BMu      = Sqr(2.0 * susyMu);
   const double scale    = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = Electroweak_constants::mtoprun * root2 / (vev * sinBeta);
   Yd(2,2) = Electroweak_constants::mbrun   * root2 / (vev * cosBeta);
   Ye(2,2) = Electroweak_constants::mtau    * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   model.set_input_parameters(input);
   model.set_scale(scale);
   model.set_loops(1);
   model.set_thresholds(3);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   ensure_tree_level_ewsb(model);
   model.calculate_DRbar_masses();
}

void setup_CMSSM_const_non_3rd_gen(CMSSM<Two_scale>& model,
                                   const CMSSM_input_parameters& input)
{
   setup_CMSSM_const(model, input);

   const double ymu = 0.1;

   Eigen::Matrix<double,3,3> Ye = model.get_Ye();
   Ye(1,1) = ymu;
   model.set_Ye(Ye);
}

void setup_data(const CMSSM<Two_scale>& model,
                Weinberg_angle::Data& data)
{
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;
   const double mt_pole = Electroweak_constants::PMTOP;
   const double scale   = model.get_scale();
   const double gY      = model.get_g1() * Sqrt(0.6);
   const double g2      = model.get_g2();
   const double e_drbar = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));

   data.fermi_contant  = Electroweak_constants::gfermi;
   data.mw_pole        = mw_pole;
   data.mz_pole        = mz_pole;
   data.mt_pole        = mt_pole;
   data.scale          = scale;
   data.gY             = gY;
   data.g2             = g2;
   data.g3             = model.get_g3();
   data.alpha_em_drbar = Sqr(e_drbar) / (4.0 * Pi);
   data.ymu            = Re(model.get_Ye(1,1));
   data.tan_beta       = model.get_vu() / model.get_vd();
   data.hmix_12        = model.get_ZH(0,1);
   data.mh_drbar       = model.get_Mhh(0);
   data.mn_drbar       = model.get_MChi();
   data.mc_drbar       = model.get_MCha();
   data.zn             = model.get_ZN();
   data.um             = model.get_UM();
   data.up             = model.get_UP();

   double msel  = 0.;
   double msmul = 0.;
   double msve  = 0.;
   double msvm  = 0.;
   const auto MSe = model.get_MSe();
   const auto ZE  = model.get_ZE();
   const auto MSv = model.get_MSv();
   const auto ZV  = model.get_ZV();

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      msel  += AbsSqr(ZE(i,0))*MSe(i);
      msmul += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      msve += AbsSqr(ZV(i,0))*MSv(i);
      msvm += AbsSqr(ZV(i,1))*MSv(i);
   }

   data.msel_drbar  = msel;
   data.msmul_drbar = msmul;
   data.msve_drbar  = msve;
   data.msvm_drbar  = msvm;

   const double pizztMZ     = Re(model.self_energy_VZ(mz_pole));
   const double piwwtMW     = Re(model.self_energy_VWm(mw_pole));
   const double piwwt0      = Re(model.self_energy_VWm(0.));
   double pizztMZ_corrected = pizztMZ;
   double piwwtMW_corrected = piwwtMW;
   double piwwt0_corrected  = piwwt0;

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = model.get_MFu(2);
   se_data.mb_drbar = model.get_MFd(2);
   se_data.gY       = gY;
   se_data.g2       = g2;

   if (model.get_thresholds() > 1) {
      pizztMZ_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole, se_data);
      piwwtMW_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(piwwtMW, mw_pole, se_data);
      piwwt0_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);
   }

   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
}

BOOST_AUTO_TEST_CASE( test_rho_2 )
{
   double r;

   r = 0.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r), CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.8;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r), CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.9;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r), CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.0;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r), CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r), CMSSM_weinberg_angle::rho_2(r), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   double delta_vb_1 = Weinberg_angle::calculate_delta_vb(outrho, outsin, data);

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = Electroweak_constants::PMTOP;
   CMSSM_weinberg_angle wein(&model, sm_parameters);
   double delta_vb_2 = wein.calculate_delta_vb(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_vb_1, delta_vb_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   double delta_r_1 = Weinberg_angle::calculate_delta_r(outrho, outsin, data);

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = Electroweak_constants::PMTOP;
   CMSSM_weinberg_angle wein(&model, sm_parameters);
   double delta_r_2 = wein.calculate_delta_r(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_r_1, delta_r_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_rho )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   double delta_rho_1 = Weinberg_angle::calculate_delta_rho(outrho, outsin, data);

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = Electroweak_constants::PMTOP;
   CMSSM_weinberg_angle wein(&model, sm_parameters);
   double delta_rho_2 = wein.calculate_delta_rho(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_rho_1, delta_rho_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_sin_theta )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   Stopwatch stopwatch;
   const double tol = 1.0e-8;
   const int maxTries = 20;
   const double rho_start = 1.0, sin_start = 0.48;

   BOOST_MESSAGE("running calculation with Weinberg_angle ...");
   stopwatch.start();
   Weinberg_angle::Data data;
   setup_data(model, data);
   Weinberg_angle wein;
   wein.set_data(data);
   wein.set_number_of_iterations(maxTries);
   wein.set_precision_goal(tol);
   int error = wein.calculate(rho_start, sin_start);
   double sin_theta_1 = wein.get_sin_theta();
   stopwatch.stop();
   double time_1 = stopwatch.get_time_in_seconds();
   BOOST_REQUIRE(error == 0);
   BOOST_MESSAGE("calculation with Weinberg_angle finished!");
   BOOST_MESSAGE("it took " << time_1 << " seconds");

   BOOST_MESSAGE("running calculation with CMSSM_weinberg_angle ...");
   stopwatch.start();
   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = Electroweak_constants::PMTOP;
   CMSSM_weinberg_angle pwein(&model, sm_parameters);
   pwein.set_number_of_iterations(maxTries);
   pwein.set_precision_goal(tol);
   double sin_theta_2;
   BOOST_REQUIRE_NO_THROW(sin_theta_2 = pwein.calculate(rho_start, sin_start));
   stopwatch.stop();
   double time_2 = stopwatch.get_time_in_seconds();
   BOOST_MESSAGE("calculation with CMSSM_weinberg_angle finished!");
   BOOST_MESSAGE("it took " << time_2 << " seconds");

   BOOST_CHECK_CLOSE_FRACTION(sin_theta_1, sin_theta_2, 1.0e-10);
}
