#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle

#include <boost/test/unit_test.hpp>

#include "conversion.hpp"
#include "softsusy.h"
#include "CMSSM_two_scale_model.hpp"
#include "test_CMSSM.hpp"

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const(fs, ss, input);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaMsbar = ss.displayDataSet().displayAlpha(ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);

   const double ss_delta_vb =
      ss.deltaVb(outrho, outsin, alphaDrbar, 0. /* ignored */);

   Weinberg_angle weinberg;

   const double gY          = ss.displayGaugeCoupling(1) * sqrt(0.6);
   const double g2          = ss.displayGaugeCoupling(2);
   const double hmu         = ss.displayYukawaElement(YE, 2, 2);
   const double mselL       = ss.displayDrBarPars().me(1, 1);
   const double msmuL       = ss.displayDrBarPars().me(1, 2);
   const double msnue       = ss.displayDrBarPars().msnu(1);
   const double msnumu      = ss.displayDrBarPars().msnu(2);
   const DoubleVector mneut = ss.displayDrBarPars().mnBpmz;
   const ComplexMatrix n    = ss.displayDrBarPars().nBpmz;
   const DoubleVector mch   = ss.displayDrBarPars().mchBpmz;
   const ComplexMatrix u    = ss.displayDrBarPars().uBpmz;
   const ComplexMatrix v    = ss.displayDrBarPars().vBpmz;

   double fs_delta_vb =
      weinberg.calculate_delta_vb(
         scale,
         outrho,
         outsin,
         mw_pole,
         mz_pole,
         alphaDrbar,
         gY,                 // displayGaugeCoupling(1) * sqrt(0.6)
         g2,                 // displayGaugeCoupling(2)
         hmu,                // = displayYukawaElement(YE, 2, 2)
         mselL,              // tree.me(1, 1)
         msmuL,              // tree.me(1, 2)
         msnue,              // tree.msnu(1)
         msnumu,             // tree.msnu(2)
         mneut,              // tree.mnBpmz
         n,                  // tree.nBpmz
         mch,                // tree.mchBpmz
         u,                  // tree.uBpmz
         v                   // tree.vBpmz
      );

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 1.0e-10);

   const auto MSe(fs.get_MSe());
   const auto ZE(fs.get_ZE());
   const auto MSv(fs.get_MSv());
   const auto ZV(fs.get_ZV());

   const double fs_gY          = fs.get_g1() * sqrt(0.6);
   const double fs_g2          = fs.get_g2();
   const double fs_hmu         = fs.get_Ye(1,1);
   double fs_mselL             = 0.;
   double fs_msmuL             = 0.;
   double fs_msnue             = 0.;
   double fs_msnumu            = 0.;
   const DoubleVector fs_mneut = ToDoubleVector(fs.get_MChi());
   const ComplexMatrix fs_n    = ToComplexMatrix(fs.get_ZN());
   const DoubleVector fs_mch   = ToDoubleVector(fs.get_MCha());
   const ComplexMatrix fs_u    = ToComplexMatrix(fs.get_UM());
   const ComplexMatrix fs_v    = ToComplexMatrix(fs.get_UP());

   for (int i = 0; i < 6; i++) {
      fs_mselL += AbsSqr(ZE(i,0))*MSe(i);
      fs_msmuL += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < 3; i++) {
      fs_msnue  += AbsSqr(ZV(i,0))*MSv(i);
      fs_msnumu += AbsSqr(ZV(i,1))*MSv(i);
   }

   BOOST_CHECK_CLOSE_FRACTION(fs_mselL , mselL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msmuL , msmuL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnue , msnue , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnumu, msnumu, 1.0e-10);

   fs_delta_vb =
      weinberg.calculate_delta_vb(
         scale,
         outrho,
         outsin,
         mw_pole,
         mz_pole,
         alphaDrbar,
         fs_gY,                 // displayGaugeCoupling(1) * sqrt(0.6)
         fs_g2,                 // displayGaugeCoupling(2)
         fs_hmu,                // = displayYukawaElement(YE, 2, 2)
         fs_mselL,              // tree.me(1, 1)
         fs_msmuL,              // tree.me(1, 2)
         fs_msnue,              // tree.msnu(1)
         fs_msnumu,             // tree.msnu(2)
         fs_mneut,              // tree.mnBpmz
         fs_n,                  // tree.nBpmz
         fs_mch,                // tree.mchBpmz
         fs_u,                  // tree.uBpmz
         fs_v                   // tree.vBpmz
      );

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 1.0e-8);
}

BOOST_AUTO_TEST_CASE( test_rho_sinTheta )
{
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const(fs, ss, input);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double tol = 1.0e-8;
   const int maxTries = 20;
   const double gfermi = Electroweak_constants::gfermi;
   const double pizztMZ = ss.piZZT(mz_pole, scale, true);
   const double piwwt0  = ss.piWWT(0., scale, true);
   const double piwwtMW = ss.piWWT(mw_pole, scale, true);
   const double alphaMsbar = ss.displayDataSet().displayAlpha(ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);

   softsusy::GMU = gfermi;

   ss.rhohat(outrho, outsin, alphaDrbar, pizztMZ, piwwt0, piwwtMW,
             tol, maxTries);

   Weinberg_angle weinberg;

   weinberg.set_number_of_iterations(maxTries);
   weinberg.set_precision_goal(tol);
   weinberg.set_alpha_em_drbar(alphaDrbar);
   weinberg.set_fermi_contant(gfermi);
   weinberg.set_self_energy_z_at_mz(pizztMZ);
   weinberg.set_self_energy_z_at_0(piwwt0);
   weinberg.set_self_energy_w_at_mw(piwwtMW);

   const double fs_sintheta = weinberg.calculate();
   const double fs_rhohat   = weinberg.get_rho_hat();

   BOOST_CHECK_CLOSE_FRACTION(outsin, fs_sintheta, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(outrho, fs_rhohat  , 1.0e-10);
}
