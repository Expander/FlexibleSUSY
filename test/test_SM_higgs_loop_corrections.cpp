
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_higgs_loop_corrections

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "pv.hpp"
#include "SM_two_scale_model.hpp"
#include "sm_twoloophiggs.hpp"

using namespace flexiblesusy;
using namespace passarino_veltman;

BOOST_AUTO_TEST_CASE( test_SM_1loop_alpha_t )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   const double mt = m.get_MFu(2);
   const double mt2 = Sqr(mt);
   const double yt = m.get_Yu(2,2);
   const double p = 125.;
   const double p2 = Sqr(p);
   const double scale = mt;
   const double scale2 = Sqr(scale);
   const double v = m.get_v();

   const double se_smh = self_energy_higgs_1loop_at_sm(
      p, scale, mt, yt);

   // top loop for p = 0, Drees p. 8
   const double se_fs_eff =
      -6 * Sqr(yt) * oneOver16PiSqr * (
         ReA0(mt2, scale2)
         + 2*mt2*ReB0(p2, mt2, mt2, scale2));

   const double se_fs =
      se_fs_eff
      + 3.*Sqr(yt)*p2*ReB0(p2,mt2,mt2,scale2) * oneOver16PiSqr;

   // tadpole
   const double t_fs =
      6 * oneOver16PiSqr * (2*yt/Sqrt(2)) * mt / v * ReA0(mt2, scale2);

   BOOST_CHECK_CLOSE_FRACTION(se_smh, se_fs + t_fs, 1.0e-10);
}
