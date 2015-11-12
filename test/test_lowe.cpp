#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowe

#include <boost/test/unit_test.hpp>
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
