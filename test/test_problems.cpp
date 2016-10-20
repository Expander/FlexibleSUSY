
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_problems

#include <boost/test/unit_test.hpp>

#include "problems.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_initialization )
{
   const std::array<std::string,3> names = {"", "", ""};
   Problems<3> problems(names);

   BOOST_CHECK(!problems.is_running_tachyon(0));
   BOOST_CHECK(!problems.is_running_tachyon(1));
   BOOST_CHECK(!problems.is_running_tachyon(2));
   BOOST_CHECK(!problems.have_tachyon());

   BOOST_CHECK(!problems.no_ewsb());
   BOOST_CHECK(!problems.no_convergence());
   BOOST_CHECK(!problems.no_perturbative());
}
