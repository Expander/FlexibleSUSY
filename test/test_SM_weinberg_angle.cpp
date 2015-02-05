#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle

#include <boost/test/unit_test.hpp>

#include "conversion.hpp"
#include "SM_two_scale_model.hpp"
#include "test_SM.hpp"
#include "logger.hpp"

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

#define ROOT2 Electroweak_constants::root2

double calculate_delta_vb_sm(
   double rho,
   double sinThetaW,
   const Weinberg_angle::Data& data
)
{
  const double g       = data.g2;
  const double gp      = data.gY;
  const double mz      = data.mz_pole;
  const double mw      = data.mw_pole;
  const double costh   = mw / mz;
  const double cw2     = Sqr(costh);
  const double sw2     = 1.0 - cw2;
  const double sinThetaW2 = Sqr(sinThetaW);
  const double outcos  = Sqrt(1.0 - sinThetaW2);
  const double q       = data.scale;
  const double alphaDRbar = data.alpha_em_drbar;
  const double hmu     = data.ymu;

  const double deltaVbSm =
     rho * alphaDRbar / (4.0 * Pi * sinThetaW2) *
     (6.0 + log(cw2) / sw2 *
      (3.5 - 2.5 * sw2 - sinThetaW2 * (5.0 - 1.5 * cw2 / Sqr(outcos))));

  return deltaVbSm;
}

double calculate_delta_r_sm(
   double rho,
   double sinThetaW,
   const Weinberg_angle::Data& data
)
{
   const double outcos = Cos(ArcSin(sinThetaW));
   const double mz = data.mz_pole;
   const double mw = data.mw_pole;
   const double mt = data.mt_pole;
   const double mh = data.mh_drbar;
   const double xt = 3.0 * data.fermi_contant * Sqr(mt) * ROOT2 * oneOver16PiSqr;
   const double alphaDRbar = data.alpha_em_drbar;
   const double g3 = data.g3;
   const double pizztMZ = data.self_energy_z_at_mz;
   const double piwwt0 = data.self_energy_w_at_0;

   const double dvb = ::calculate_delta_vb_sm(rho, sinThetaW, data);

   const double deltaR = rho * piwwt0 / Sqr(mw) -
      pizztMZ / Sqr(mz) + dvb;

   const double deltaR2LoopSm = alphaDRbar * Sqr(g3) /
      (16.0 * Sqr(Pi) * Pi * Sqr(sinThetaW) * Sqr(outcos)) *
      (2.145 * Sqr(mt) / Sqr(mz) + 0.575 * log(mt / mz) - 0.224
       - 0.144 * Sqr(mz) / Sqr(mt)) -
      Sqr(xt) *
      Weinberg_angle::rho_2(mh / mt) * (1.0 - deltaR) * rho / 3.0;

   const double deltaR_full = deltaR + deltaR2LoopSm;

   return deltaR_full;
}

double calculate_delta_rho_sm(
   double rho,
   double sinThetaW,
   const Weinberg_angle::Data& data
)
{
   const double mz = data.mz_pole;
   const double mw = data.mw_pole;
   const double mt = data.mt_pole;
   const double mh = data.mh_drbar;
   const double xt = 3.0 * data.fermi_contant * Sqr(mt) * ROOT2 * oneOver16PiSqr;
   const double alphaDRbar = data.alpha_em_drbar;
   const double g3 = data.g3;
   const double pizztMZ = data.self_energy_z_at_mz;
   const double piwwtMW = data.self_energy_w_at_mw;

   const double deltaRho2LoopSm = alphaDRbar * Sqr(g3) /
      (16.0 * Pi * Sqr(Pi) * Sqr(sinThetaW)) *
      (-2.145 * Sqr(mt) / Sqr(mw) + 1.262 * log(mt / mz) - 2.24
       - 0.85 * Sqr(mz)
       / Sqr(mt)) + Sqr(xt) *
      Weinberg_angle::rho_2(mh / mt) / 3.0;

   const double deltaRhoOneLoop = pizztMZ / (rho * Sqr(mz))
      - piwwtMW / Sqr(mw);

   const double deltaRho = deltaRhoOneLoop + deltaRho2LoopSm;

   return deltaRho;
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   SM<Two_scale> fs;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(fs, input);

   const double scale   = fs.get_scale();
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaDrbar = Electroweak_constants::aem;

   const double fs_gY  = fs.get_g1() * sqrt(0.6);
   const double fs_g2  = fs.get_g2();
   const double fs_hmu = fs.get_Ye(1,1);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.gY = fs_gY;
   data.g2 = fs_g2;
   data.ymu = fs_hmu;

   const double ss_delta_vb =
      ::calculate_delta_vb_sm(outrho, outsin, data);

   const double fs_delta_vb =
      Weinberg_angle::calculate_delta_vb_sm(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   SM<Two_scale> fs;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(fs, input);

   const double scale   = fs.get_scale();
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaDrbar = Electroweak_constants::aem;

   const double fs_gY  = fs.get_g1() * sqrt(0.6);
   const double fs_g2  = fs.get_g2();
   const double fs_hmu = fs.get_Ye(1,1);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = Electroweak_constants::gfermi;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = Electroweak_constants::PMTOP;
   data.mh_drbar = Electroweak_constants::MH;
   data.gY = fs_gY;
   data.g2 = fs_g2;
   data.ymu = fs_hmu;

   const double ss_delta_r =
      ::calculate_delta_r_sm(outrho, outsin, data);

   const double fs_delta_r =
      Weinberg_angle::calculate_delta_r(outrho, outsin, data, false);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_r, fs_delta_r, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_rho )
{
   SM<Two_scale> fs;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(fs, input);

   const double scale   = fs.get_scale();
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaDrbar = Electroweak_constants::aem;

   const double fs_gY  = fs.get_g1() * sqrt(0.6);
   const double fs_g2  = fs.get_g2();
   const double fs_hmu = fs.get_Ye(1,1);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = Electroweak_constants::gfermi;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = Electroweak_constants::PMTOP;
   data.mh_drbar = Electroweak_constants::MH;
   data.gY = fs_gY;
   data.g2 = fs_g2;
   data.ymu = fs_hmu;

   const double ss_delta_r =
      ::calculate_delta_rho_sm(outrho, outsin, data);

   const double fs_delta_r =
      Weinberg_angle::calculate_delta_rho(outrho, outsin, data, false);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_r, fs_delta_r, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_rho_sinTheta )
{
   SM<Two_scale> fs;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(fs, input);

   const double scale   = fs.get_scale();
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaDrbar = Electroweak_constants::aem;

   const double fs_gY  = fs.get_g1() * sqrt(0.6);
   const double fs_g2  = fs.get_g2();
   const double fs_hmu = fs.get_Ye(1,1);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = Electroweak_constants::gfermi;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = Electroweak_constants::PMTOP;
   data.mh_drbar = Electroweak_constants::MH;
   data.gY = fs_gY;
   data.g2 = fs_g2;
   data.ymu = fs_hmu;

   Weinberg_angle weinberg;
   weinberg.disable_susy_contributions();
   weinberg.set_number_of_iterations(20);
   weinberg.set_precision_goal(1.0e-8);
   weinberg.set_data(data);

   const double fs_sin = weinberg.calculate();

   BOOST_MESSAGE("SM sin(ThetaW) = " << fs_sin);
}
