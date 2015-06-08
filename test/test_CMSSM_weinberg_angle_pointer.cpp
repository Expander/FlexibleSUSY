#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle_pointer

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

#define private public

#include "weinberg_angle.hpp"

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

   double mselL  = 0.;
   double msmuL  = 0.;
   double msnue  = 0.;
   double msnumu = 0.;
   const auto MSe = model.get_MSe();
   const auto ZE  = model.get_ZE();
   const auto MSv = model.get_MSv();
   const auto ZV  = model.get_ZV();

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      mselL += AbsSqr(ZE(i,0))*MSe(i);
      msmuL += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      msnue  += AbsSqr(ZV(i,0))*MSv(i);
      msnumu += AbsSqr(ZV(i,1))*MSv(i);
   }

   data.msel_drbar  = mselL;
   data.msmul_drbar = msmuL;
   data.msve_drbar  = msnue;
   data.msvm_drbar  = msnumu;

   const double pizztMZ     = Re(model.self_energy_VZ(mz_pole));
   const double piwwtMW     = Re(model.self_energy_VWm(mw_pole));
   const double piwwt0      = Re(model.self_energy_VWm(0.));
   double pizztMZ_corrected = 0.;
   double piwwtMW_corrected = 0.;
   double piwwt0_corrected  = 0.;

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
   BOOST_REQUIRE(1 > 0);
/*   Weinberg_angle wein;
   double r;

   r = 0.1;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 1.8;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 1.9;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 2.0;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 2.1;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);*/
}

/*
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
   double fs_delta_vb =
      Weinberg_angle::calculate_delta_vb(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 3.0e-9);
}
*/

/*
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
   double fs_delta_r =
      Weinberg_angle::calculate_delta_r(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_r, fs_delta_r, 3.0e-9);
}
*/

/*
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
   double fs_delta_rho =
      Weinberg_angle::calculate_delta_rho(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_rho, fs_delta_rho, 1.0e-10);
}
*/

/*
BOOST_AUTO_TEST_CASE( test_rho_sinTheta )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   Weinberg_angle::Data data;
   setup_data(model, data);

   const double tol = 1.0e-8;
   const int maxTries = 20;

   const double rho_start = 1.0, sin_start = 0.48;

   Weinberg_angle weinberg;
   weinberg.set_number_of_iterations(maxTries);
   weinberg.set_precision_goal(tol);
   weinberg.set_data(data);

   BOOST_MESSAGE("running Weinberg_angle::calculate() ...");

   const int error = weinberg.calculate(rho_start, sin_start);

   BOOST_REQUIRE(error == 0);
   const double fs_sintheta = weinberg.get_sin_theta();

   BOOST_MESSAGE("running Weinberg_angle::calculate() finished!");

   const double fs_rhohat = weinberg.get_rho_hat();

   BOOST_CHECK_CLOSE_FRACTION(outsin, fs_sintheta, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(outrho, fs_rhohat  , 1.0e-10);
}
*/

/*
BOOST_AUTO_TEST_CASE( test_self_energy_top_correction )
{
   CMSSM<Two_scale> model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   Weinberg_angle::Data data;
   setup_data(model, data);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();
   const double mt_drbar = fs.get_MFu(2);
   const double mb_drbar = fs.get_MFd(2);

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);
   BOOST_CHECK_EQUAL(fs.get_scale(), mz_pole);
   BOOST_CHECK_CLOSE_FRACTION(mt_drbar, ss.displayDrBarPars().mt, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(mb_drbar, ss.displayDrBarPars().mb, 1.0e-10);

   // reference values
   const double pizztMZ = ss.piZZT(mz_pole, scale, true);
   const double piwwtMW = ss.piWWT(mw_pole, scale, true);
   const double piwwt0  = ss.piWWT(0.     , scale, true);

   const double ss_pizztMZ = ss.piZZT(mz_pole, scale, false);
   const double ss_piwwtMW = ss.piWWT(mw_pole, scale, false);
   const double ss_piwwt0  = ss.piWWT(0.     , scale, false);
   const double fs_pizztMZ = Re(fs.self_energy_VZ(mz_pole));
   const double fs_piwwtMW = Re(fs.self_energy_VWm(mw_pole));
   const double fs_piwwt0  = Re(fs.self_energy_VWm(0));

   BOOST_CHECK_CLOSE_FRACTION(ss_pizztMZ, fs_pizztMZ, 5.0e-08);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwtMW, fs_piwwtMW, 3.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwt0 , fs_piwwt0 , 5.0e-04);

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale = scale;
   se_data.mt_pole = data.mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY = data.gY;
   se_data.g2 = data.g2;

   const double fs_pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(fs_pizztMZ, mz_pole, se_data);

   const double fs_piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(fs_piwwtMW, mw_pole, se_data);

   const double fs_piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(fs_piwwt0, 0., se_data);

   BOOST_CHECK_CLOSE_FRACTION(fs_pizztMZ_corrected, pizztMZ, 5.0e-08);
   BOOST_CHECK_CLOSE_FRACTION(fs_piwwtMW_corrected, piwwtMW, 3.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(fs_piwwt0_corrected , piwwt0 , 4.0e-04);
}
*/
