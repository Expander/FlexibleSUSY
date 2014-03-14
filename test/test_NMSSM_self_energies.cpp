
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_self_energies

#include <boost/test/unit_test.hpp>

#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;


BOOST_AUTO_TEST_CASE( test_NMSSM_self_energy_CP_even_higgs )
{
   NMSSM_input_parameters input;
   input.m0 = 300.; // avoids tree-level tachyons
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   // initial guess
   m.set_Kappa(0.1);
   m.set_vS(5000.);
   m.set_ms2(-Sqr(input.m0));
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   s.setKappa(m.get_Kappa());
   s.setSvev(m.get_vS());
   s.setMsSquared(m.get_ms2());
   s.setMh1Squared(m.get_mHd2());
   s.setMh2Squared(m.get_mHu2());

   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   // check tree-level
   DoubleVector hh_ss(s.displayDrBarPars().mh0);
   Eigen::Array<double,3,1> hh_fs(m.get_Mhh());

   BOOST_CHECK_CLOSE_FRACTION(hh_ss(1), hh_fs(0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(2), hh_fs(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(3), hh_fs(2), 1.0e-10);

   // solve EWSB
   int loop_order = 2;
   softsusy::numRewsbLoops = loop_order;
   softsusy::numHiggsMassLoops = loop_order;
   m.set_ewsb_loop_order(loop_order);
   m.set_pole_mass_loop_order(loop_order);

   ensure_n_loop_ewsb(m, loop_order);
   ensure_n_loop_ewsb(s, loop_order);

   const double kappa_ss = s.displayKappa();
   const double vS_ss    = s.displaySvev();
   const double ms2_ss   = s.displayMsSquared();

   const double kappa_fs = m.get_Kappa();
   const double vS_fs    = m.get_vS();
   const double ms2_fs   = m.get_ms2();

   BOOST_CHECK_CLOSE_FRACTION(kappa_ss, kappa_fs, 1.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(vS_ss   , vS_fs   , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(ms2_ss  , ms2_fs  , 1.0e-10);

   const int accuracy = loop_order;
   const double piWWT = 0., pizztMS = 0.;

   // setting initial conditions for the iteration
   softsusy::sPhysical physical(s.displayDrBarPars());
   physical.mh0(1) = s.displayDrBarPars().mh0(1);
   physical.mh0(2) = s.displayDrBarPars().mh0(2);
   physical.mh0(3) = s.displayDrBarPars().mh0(3);
   physical.mA0(1) = s.displayDrBarPars().mA0(1);
   physical.mA0(2) = s.displayDrBarPars().mA0(2);

   s.setPhys(physical); // initialization
   s.higgs(accuracy, piWWT, pizztMS, physical); // does one iteration
   s.setPhys(physical); // now modified

   m.set_number_of_mass_iterations(1);
   m.calculate_Mhh_pole();

   hh_ss = s.displayPhys().mh0;
   hh_fs = m.get_physical().Mhh;

   BOOST_CHECK_CLOSE_FRACTION(hh_ss(1), hh_fs(0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(2), hh_fs(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(3), hh_fs(2), 1.0e-10);
}
