#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloop_mt

#include <boost/test/unit_test.hpp>
#include "mssm_twoloop_mt.hpp"
#include "wrappers.hpp"
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::mssm_twoloop_mt;

const double g3 = std::sqrt(4*Pi* 0.1184);
const double mt = 173.34;

BOOST_AUTO_TEST_CASE(test_universal)
{
   std::vector<Parameters> pars = {
      // g3 mt mg mst1 mst2 msusy Xt Q
      { g3, mt, 1000., 1000., 1000., 1000.,                  0., 1000. },
      { g3, mt, 1000., 1000., 1000., 1000.,                  0., 1100. },
      { g3, mt, 1000., 1000., 1000., 1000.,  std::sqrt(6.)*1000, 1000. },
      { g3, mt, 1000., 1000., 1000., 1000.,  std::sqrt(6.)*1000, 1100. },
      { g3, mt, 1000., 1000., 1000., 1000., -std::sqrt(6.)*1000, 1000. },
      { g3, mt, 1000., 1000., 1000., 1000., -std::sqrt(6.)*1000, 1100. }
   };

   std::cout << "# universal SUSY masses" << std::endl;
   for (const auto p: pars)
      std::cout << dMt_over_mt_2loop(p) << std::endl;
}

BOOST_AUTO_TEST_CASE(test_non_universal)
{
   std::vector<Parameters> pars = {
      // g3 mt mg mst1 mst2 msusy Xt Q
      { g3, mt, 1000., 1200., 1300., 1400.,                  0., 1000. },
      { g3, mt, 1000., 1200., 1300., 1400.,                  0., 1100. },
      { g3, mt, 1000., 1200., 1300., 1400.,  std::sqrt(6.)*1300, 1000. },
      { g3, mt, 1000., 1200., 1300., 1400.,  std::sqrt(6.)*1300, 1100. },
      { g3, mt, 1000., 1200., 1300., 1400., -std::sqrt(6.)*1300, 1000. },
      { g3, mt, 1000., 1200., 1300., 1400., -std::sqrt(6.)*1300, 1100. }
   };

   std::cout << "# non-universal SUSY masses" << std::endl;
   for (const auto p: pars)
      std::cout << dMt_over_mt_1loop_qcd(p) << ' '
                << dMt_over_mt_1loop(p) << ' '
                << dMt_over_mt_2loop_qcd(p) << ' '
                << dMt_over_mt_2loop(p) << std::endl;
}
