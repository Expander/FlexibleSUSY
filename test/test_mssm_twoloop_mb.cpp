#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloop_mb

#include <boost/test/unit_test.hpp>
#include "mssm_twoloop_mb.hpp"
#include "wrappers.hpp"
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::mssm_twoloop_mb;

const double g3 = std::sqrt(4*Pi* 0.1184);
const double mt = 173.34;
const double mb = 3.20;

BOOST_AUTO_TEST_CASE(test_non_universal)
{
   std::vector<Parameters> pars = {
      // g3 mt  mb     mg   mst1   mst2   msb1   msb2   msusy                   Xt  Xb      Q
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,                  0., 0., 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,                  0., 0., 1100. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,  std::sqrt(6.)*1300, 0., 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,  std::sqrt(6.)*1300, 0., 1100. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400., -std::sqrt(6.)*1300, 0., 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400., -std::sqrt(6.)*1300, 0., 1100. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,                  0.,  std::sqrt(6.), 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,                  0., -std::sqrt(6.), 1100. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,  std::sqrt(6.)*1300,  std::sqrt(6.), 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400.,  std::sqrt(6.)*1300, -std::sqrt(6.), 1100. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400., -std::sqrt(6.)*1300,  std::sqrt(6.), 1000. },
      { g3, mt, mb, 1000., 1200., 1300., 1500., 1600.,  1400., -std::sqrt(6.)*1300, -std::sqrt(6.), 1100. }
   };

   std::cout << "# non-universal SUSY masses" << std::endl;
   for (const auto p: pars)
      std::cout << delta_mb_2loop(p) << std::endl;
}
