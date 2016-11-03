#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowe

#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"
#include "lowe.h"
#include "conversion.hpp"

using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_to )
{
   QedQcd lowe_Mz;
   lowe_Mz.setPoleMt(173.5);
   lowe_Mz.setAlpha(ALPHAS, 0.118); // running alpha_s at current scale
   lowe_Mz.setAlphaSInput(0.118);   // input alpha_s(MZ)
   QedQcd lowe_Mz_new(lowe_Mz);

   lowe_Mz.toMz(); // uses running alpha_s at current scale
   lowe_Mz_new.to(lowe_Mz.displayPoleMZ()); // uses input alpha_s(MZ)

   BOOST_CHECK_CLOSE(lowe_Mz.displayMbMb(), lowe_Mz_new.displayMbMb(), 1e-10);
   BOOST_CHECK_CLOSE(lowe_Mz.displayAlpha(ALPHA) , lowe_Mz_new.displayAlpha(ALPHA) , 1e-10);
   BOOST_CHECK_CLOSE(lowe_Mz.displayAlpha(ALPHAS), lowe_Mz_new.displayAlpha(ALPHAS), 1e-10);

   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mUp)      , lowe_Mz_new.displayMass(mUp)      , 1.1e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mDown)    , lowe_Mz_new.displayMass(mDown)    , 1.1e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mStrange) , lowe_Mz_new.displayMass(mStrange) , 1.1e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mCharm)   , lowe_Mz_new.displayMass(mCharm)   , 5e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mBottom)  , lowe_Mz_new.displayMass(mBottom)  , 0.5);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mElectron), lowe_Mz_new.displayMass(mElectron), 0.5);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mMuon)    , lowe_Mz_new.displayMass(mMuon)    , 0.5);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mTau)     , lowe_Mz_new.displayMass(mTau)     , 0.5);

   BOOST_CHECK_CLOSE(lowe_Mz.displayPoleMb(), lowe_Mz_new.displayPoleMb(), 0.5);

   BOOST_TEST_MESSAGE(lowe_Mz);
   BOOST_TEST_MESSAGE(lowe_Mz_new);
}

BOOST_AUTO_TEST_CASE( test_to_recall )
{
   QedQcd lowe_Mz;
   lowe_Mz.setPoleMt(173.5);
   lowe_Mz.setAlpha(ALPHAS, 0.118);
   lowe_Mz.setMu(lowe_Mz.displayPoleMZ());
   lowe_Mz.to(100.);

   QedQcd lowe_Mz_new(lowe_Mz);
   lowe_Mz.setMu(100.);
   lowe_Mz_new.to(100.);

   BOOST_CHECK_LT(flexiblesusy::MaxRelDiff(lowe_Mz.display_input(), lowe_Mz_new.display_input()), 1e-10);
   BOOST_CHECK_LT(flexiblesusy::MaxRelDiff(flexiblesusy::ToEigenArray(lowe_Mz.display()),
                                           flexiblesusy::ToEigenArray(lowe_Mz_new.display())), 1e-10);

   BOOST_TEST_MESSAGE(lowe_Mz);
   BOOST_TEST_MESSAGE(lowe_Mz_new);
}
