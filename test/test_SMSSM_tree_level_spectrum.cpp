
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_SMSSM_tree_level_masses )
{
   SMSSM_input_parameters input;
   input.m0 = 250.;
   SMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_SMSSM(m, s, input);

   m.calculate_DRbar_parameters();
   s.calcDrBarPars();

   // neutral CP even Higgs
   const DoubleVector hh(ToDoubleVector(m.get_Mhh()));
   const DoubleVector mh0(s.displayDrBarPars().mh0);
   BOOST_CHECK_CLOSE(hh(1), mh0(1), 1.0e-12);
   BOOST_CHECK_CLOSE(hh(2), mh0(2), 1.0e-12);
   BOOST_CHECK_CLOSE(hh(3), mh0(3), 1.0e-12);

   // neutral CP odd Higgs
   const DoubleVector Ah(ToDoubleVector(m.get_MAh()));
   const double mz = s.displayMzRun();
   DoubleVector mA0(3);
   mA0(1) = std::min(s.displayDrBarPars().mA0(1), mz);
   mA0(2) = std::max(s.displayDrBarPars().mA0(1), mz);
   mA0(3) = s.displayDrBarPars().mA0(2);
   BOOST_CHECK_CLOSE(Ah(1), mA0(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Ah(2), mA0(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Ah(3), mA0(3), 1.0e-12);

   // charged Higgs
   const DoubleVector Hpm(ToDoubleVector(m.get_MHpm()));
   const double mHpm = s.displayDrBarPars().mHpm;
   const double mw = s.displayMwRun();
   BOOST_CHECK_CLOSE(Hpm(1), mw  , 1.0e-12);
   BOOST_CHECK_CLOSE(Hpm(2), mHpm, 1.0e-12);

   // neutralinos
   const DoubleVector mneut(ToDoubleVector(m.get_MChi()));
   const DoubleVector mnBpmz = s.displayDrBarPars().mnBpmz;
   BOOST_CHECK_CLOSE(mneut(1), mnBpmz(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(2), mnBpmz(2), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(3), mnBpmz(3), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(4), mnBpmz(4), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(5), mnBpmz(5), 1.0e-12);

   // charginos
   const DoubleVector mch(ToDoubleVector(m.get_MCha()));
   const DoubleVector mchBpmz = s.displayDrBarPars().mchBpmz;
   BOOST_CHECK_CLOSE(mch(1), mchBpmz(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mch(2), mchBpmz(2), 1.0e-12);

   // photon, W, Z, gluon mass
   const double vp = m.get_MVP();
   const double vz = m.get_MVZ();
   const double vw = m.get_MVWm();
   const double vg = m.get_MVG();
   BOOST_CHECK_EQUAL(vp, 0.0);
   BOOST_CHECK_CLOSE(vz, mz, 1.0e-12);
   BOOST_CHECK_CLOSE(vw, mw, 1.0e-12);
   BOOST_CHECK_EQUAL(vg, 0.0);

   // down-type squarks
   const DoubleVector Sd(ToDoubleVector(m.get_MSd()));
   const DoubleVector md(s.displayDrBarPars().md.flatten().sort());
   BOOST_CHECK_CLOSE(Sd(1), md(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(2), md(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(3), md(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(4), md(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(5), md(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(6), md(6), 1.0e-12);

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(m.get_MSu()));
   const DoubleVector mu(s.displayDrBarPars().mu.flatten().sort());
   BOOST_CHECK_CLOSE(Su(1), mu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(2), mu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(3), mu(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(4), mu(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(5), mu(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(6), mu(6), 1.0e-12);

   // down-type sleptons
   const DoubleVector Se(ToDoubleVector(m.get_MSe()));
   const DoubleVector me(s.displayDrBarPars().me.flatten().sort());
   BOOST_CHECK_CLOSE(Se(1), me(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(2), me(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(3), me(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(4), me(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(5), me(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(6), me(6), 1.0e-12);

   // gluinos
   const double MGlu = m.get_MGlu();
   const double mGluino = s.displayDrBarPars().mGluino;
   BOOST_CHECK_CLOSE(MGlu, mGluino, 1.0e-12);

   // neutrinos
   const DoubleVector MFv(ToDoubleVector(m.get_MFv()));
   BOOST_CHECK_EQUAL(MFv(1), 0.0);
   BOOST_CHECK_EQUAL(MFv(2), 0.0);
   BOOST_CHECK_EQUAL(MFv(3), 0.0);

   // leptons
   const DoubleVector MFe(ToDoubleVector(m.get_MFe()));
   BOOST_CHECK_EQUAL(MFe(1), 0.0);
   BOOST_CHECK_EQUAL(MFe(2), 0.0);
   BOOST_CHECK_CLOSE(MFe(3), s.displayDrBarPars().mtau, 1.0e-12);

   // ups
   const DoubleVector MFu(ToDoubleVector(m.get_MFu()));
   BOOST_CHECK_EQUAL(MFu(1), 0.0);
   BOOST_CHECK_EQUAL(MFu(2), 0.0);
   BOOST_CHECK_CLOSE(MFu(3), s.displayDrBarPars().mt, 1.0e-12);

   // downs
   const DoubleVector MFd(ToDoubleVector(m.get_MFd()));
   BOOST_CHECK_EQUAL(MFd(1), 0.0);
   BOOST_CHECK_EQUAL(MFd(2), 0.0);
   BOOST_CHECK_CLOSE(MFd(3), s.displayDrBarPars().mb, 1.0e-12);
}
