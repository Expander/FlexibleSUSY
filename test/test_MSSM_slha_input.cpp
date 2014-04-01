
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_slha_input

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_model.hpp"
#include "MSSM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MSSM_slha_reading )
{
   MSSM<Two_scale> model;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * PI * alpha1);
   const double g2 = sqrt(4 * PI * alpha2);
   const double g3 = sqrt(4 * PI * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);

   model.set_scale(91);
   model.set_loops(2);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(sqr(m0));
   model.set_mHu2(sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   MSSM_slha_io slha_io;
   slha_io.set_spectrum(model);

   std::string slha_file("test/test_MSSM_slha_input.out.spc");
   slha_io.write_to_file(slha_file);

   model.clear();
   slha_io.clear();

   BOOST_CHECK_EQUAL(model.get_g1(), 0.);
   BOOST_CHECK_EQUAL(model.get_g2(), 0.);
   BOOST_CHECK_EQUAL(model.get_g3(), 0.);

   slha_io.read_from_file(slha_file);
   slha_io.fill(model);

   BOOST_CHECK_EQUAL(model.get_g1(), g1);
   BOOST_CHECK_EQUAL(model.get_g2(), g2);
   BOOST_CHECK_EQUAL(model.get_g3(), g3);
}
