#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle

#include <boost/test/unit_test.hpp>

#include "weinberg_angle.hpp"
#include "softsusy.h"
#include "CMSSM_two_scale_model.hpp"
#include "test_CMSSM.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

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
