
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_NMSSM_linking

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_NMSSM_linking )
{
   CMSSM<Two_scale> mssm;
   NMSSM<Two_scale> nmssm;

   BOOST_CHECK(true);
}
