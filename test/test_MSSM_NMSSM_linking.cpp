
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_NMSSM_linking

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MSSM_NMSSM_linking )
{
   MSSM<Two_scale> mssm;
   NMSSM<Two_scale> nmssm;

   BOOST_CHECK(true);
}
