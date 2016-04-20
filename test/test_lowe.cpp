#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowe

#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"
#include "lowe.h"

using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_toMz_toQ )
{
   QedQcd lowe_MZ, lowe_Q;
   lowe_MZ.toMz();
   lowe_Q.to(lowe_MZ.displayPoleMZ());

   BOOST_CHECK(lowe_MZ == lowe_Q);
}

BOOST_AUTO_TEST_CASE( test_toMt_toQ )
{
   QedQcd lowe_Mt, lowe_Q;
   lowe_Mt.toMt();
   lowe_Q.to(lowe_Mt.displayPoleMt());

   BOOST_MESSAGE(lowe_Mt);
   BOOST_MESSAGE(lowe_Q);

   BOOST_CHECK(lowe_Mt == lowe_Q);
}

BOOST_AUTO_TEST_CASE( test_toMz_recall )
{
   QedQcd lowe_Mz;
   lowe_Mz.setPoleMt(173.5);
   lowe_Mz.setAlpha(ALPHAS, 0.118);
   QedQcd lowe_Mz_new(lowe_Mz);

   lowe_Mz.toMz();
   lowe_Mz_new.to2(lowe_Mz.displayPoleMZ());

   BOOST_CHECK_CLOSE(lowe_Mz.displayMbMb(), lowe_Mz_new.displayMbMb(), 1e-10);
   BOOST_CHECK_CLOSE(lowe_Mz.displayAlpha(ALPHA) , lowe_Mz_new.displayAlpha(ALPHA) , 1e-10);
   BOOST_CHECK_CLOSE(lowe_Mz.displayAlpha(ALPHAS), lowe_Mz_new.displayAlpha(ALPHAS), 1e-10);

   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mUp)      , lowe_Mz_new.displayMass(mUp)      , 5e-3);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mDown)    , lowe_Mz_new.displayMass(mDown)    , 5e-3);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mStrange) , lowe_Mz_new.displayMass(mStrange) , 5e-3);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mCharm)   , lowe_Mz_new.displayMass(mCharm)   , 5e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mBottom)  , lowe_Mz_new.displayMass(mBottom)  , 3e-2);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mElectron), lowe_Mz_new.displayMass(mElectron), 2.);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mMuon)    , lowe_Mz_new.displayMass(mMuon)    , 2.);
   BOOST_CHECK_CLOSE(lowe_Mz.displayMass(mTau)     , lowe_Mz_new.displayMass(mTau)     , 2.);

   BOOST_CHECK_CLOSE(lowe_Mz.displayPoleMb(), lowe_Mz_new.displayPoleMb(), 5e-2);

   BOOST_MESSAGE(lowe_Mz);
   BOOST_MESSAGE(lowe_Mz_new);
}
