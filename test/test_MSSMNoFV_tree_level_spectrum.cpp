
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMNoFV_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test.h"
#include "test_MSSMNoFV.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSMNoFV_two_scale_model.hpp"

#define COMPARE_MASS(p) TEST_EQUALITY(a.get_##p(), b.get_##p());

using namespace flexiblesusy;

void test_spectrum_equality(const MSSMNoFV<Two_scale>& a, const MSSM<Two_scale>& b)
{
   COMPARE_MASS(Mhh);
   COMPARE_MASS(MAh);
   COMPARE_MASS(MHpm);

   COMPARE_MASS(MChi);
   COMPARE_MASS(MCha);
   COMPARE_MASS(MGlu);

   COMPARE_MASS(MVP);
   COMPARE_MASS(MVZ);
   COMPARE_MASS(MVWm);
   COMPARE_MASS(MVG);

   // up-type quarks
   const Eigen::Array<double,2,1> MSu(a.get_MSu());
   const Eigen::Array<double,2,1> MSc(a.get_MSc());
   const Eigen::Array<double,2,1> MSt(a.get_MSt());
   const Eigen::Array<double,6,1> MSu_full(b.get_MSu());
   Eigen::Array<double,6,1> MSu_mixed;
   MSu_mixed(0) = MSu(0);
   MSu_mixed(1) = MSu(1);
   MSu_mixed(2) = MSc(0);
   MSu_mixed(3) = MSc(1);
   MSu_mixed(4) = MSt(0);
   MSu_mixed(5) = MSt(1);
   std::sort(MSu_mixed.data(), MSu_mixed.data() + MSu_mixed.size());

   TEST_EQUALITY(MSu_full, MSu_mixed);

   // down-type quarks
   const Eigen::Array<double,2,1> MSd(a.get_MSd());
   const Eigen::Array<double,2,1> MSs(a.get_MSs());
   const Eigen::Array<double,2,1> MSb(a.get_MSb());
   const Eigen::Array<double,6,1> MSd_full(b.get_MSd());
   Eigen::Array<double,6,1> MSd_mixed;
   MSd_mixed(0) = MSd(0);
   MSd_mixed(1) = MSd(1);
   MSd_mixed(2) = MSs(0);
   MSd_mixed(3) = MSs(1);
   MSd_mixed(4) = MSb(0);
   MSd_mixed(5) = MSb(1);
   std::sort(MSd_mixed.data(), MSd_mixed.data() + MSd_mixed.size());

   TEST_EQUALITY(MSd_full, MSd_mixed);

   // up-type leptons
   const double MSveL(a.get_MSveL());
   const double MSvmL(a.get_MSvmL());
   const double MSvtL(a.get_MSvtL());
   const Eigen::Array<double,3,1> MSv_full(b.get_MSv());
   Eigen::Array<double,3,1> MSv_mixed;
   MSv_mixed(0) = MSveL;
   MSv_mixed(1) = MSvmL;
   MSv_mixed(2) = MSvtL;
   std::sort(MSv_mixed.data(), MSv_mixed.data() + MSv_mixed.size());

   TEST_EQUALITY(MSv_full, MSv_mixed);

   // down-type leptons
   const Eigen::Array<double,2,1> MSe(a.get_MSe());
   const Eigen::Array<double,2,1> MSm(a.get_MSm());
   const Eigen::Array<double,2,1> MStau(a.get_MStau());
   const Eigen::Array<double,6,1> MSe_full(b.get_MSe());
   Eigen::Array<double,6,1> MSe_mixed;
   MSe_mixed(0) = MSe(0);
   MSe_mixed(1) = MSe(1);
   MSe_mixed(2) = MSm(0);
   MSe_mixed(3) = MSm(1);
   MSe_mixed(4) = MStau(0);
   MSe_mixed(5) = MStau(1);
   std::sort(MSe_mixed.data(), MSe_mixed.data() + MSe_mixed.size());

   TEST_EQUALITY(MSe_full, MSe_mixed);
}

BOOST_AUTO_TEST_CASE( test_MSSMNoFV_tree_level_spectrum )
{
   MSSMNoFV_input_parameters input;
   MSSMNoFV<Two_scale> m1;
   MSSM<Two_scale> m2;
   setup_MSSM_models(m1, m2, input);

   test_spectrum_equality(m1, m2);
   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("spectra are not equal");
      gErrors = 0;
   }
}
