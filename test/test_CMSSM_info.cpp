
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_info

#include <boost/test/unit_test.hpp>
#include <cstdlib>

#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_info.hpp"

using namespace flexiblesusy;

double get_random() {
   return 1.0 * rand() / RAND_MAX;
}

#define TEST_SCALAR(p)                                                  \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                             \
      BOOST_CHECK_EQUAL(mssm.get_##p(), 0.0);                           \
      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());       \
      const double value = get_random();                                \
      parameters(CMSSM_info::p) = value;                                 \
      mssm.set(parameters);                                             \
      BOOST_CHECK_EQUAL(mssm.get_##p(), value);                         \
   }

#define TEST_DOUBLE_MATRIX_3x3(p)                                       \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                             \
      for (int i = 0; i < 3; i++)                                       \
         for (int k = 0; k < 3; k++)                                    \
            BOOST_CHECK_EQUAL(mssm.get_##p(i,k), 0.0);                  \
      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());       \
      Eigen::Matrix<double,3,3> p;                                      \
      p << get_random(), get_random(), get_random(),                    \
         get_random(), get_random(), get_random(),                      \
         get_random(), get_random(), get_random();                      \
      parameters(CMSSM_info::p##0_0) = p(0,0);                          \
      parameters(CMSSM_info::p##0_1) = p(0,1);                          \
      parameters(CMSSM_info::p##0_2) = p(0,2);                          \
      parameters(CMSSM_info::p##1_0) = p(1,0);                          \
      parameters(CMSSM_info::p##1_1) = p(1,1);                          \
      parameters(CMSSM_info::p##1_2) = p(1,2);                          \
      parameters(CMSSM_info::p##2_0) = p(2,0);                          \
      parameters(CMSSM_info::p##2_1) = p(2,1);                          \
      parameters(CMSSM_info::p##2_2) = p(2,2);                          \
      mssm.set(parameters);                                             \
      for (int i = 0; i < 3; i++)                                       \
         for (int k = 0; k < 3; k++)                                    \
            BOOST_CHECK_EQUAL(mssm.get_##p(i,k), p(i,k));               \
   }                                                                    \

BOOST_AUTO_TEST_CASE( test_CMSSM_parameter_enum )
{
   // This test checks that the order of the enum entries
   // CMSSM_info::Parameters matches the order of the parameters in the
   // model class (as used in the set() and get() functions).
   //
   // If this test works, then the user can fill a parameter vector
   // using the enum entries and then set all model parameters at once
   // via the set() function.

   srand(1);

   {
      CMSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(CMSSM_info::NUMBER_OF_PARAMETERS,
                        mssm.get_number_of_parameters());
   }

   // ==================== susy parameters ====================

   TEST_SCALAR(g1);
   TEST_SCALAR(g2);
   TEST_SCALAR(g3);
   TEST_DOUBLE_MATRIX_3x3(Yu);
   TEST_DOUBLE_MATRIX_3x3(Yd);
   TEST_DOUBLE_MATRIX_3x3(Ye);
   TEST_SCALAR(Mu);
   TEST_SCALAR(vd);
   TEST_SCALAR(vu);

   // ==================== soft parameters ====================

   TEST_DOUBLE_MATRIX_3x3(TYu);
   TEST_DOUBLE_MATRIX_3x3(TYd);
   TEST_DOUBLE_MATRIX_3x3(TYe);
   TEST_SCALAR(BMu);
   TEST_SCALAR(MassB);
   TEST_SCALAR(MassWB);
   TEST_SCALAR(MassG);
   TEST_SCALAR(mHd2);
   TEST_SCALAR(mHu2);
   TEST_DOUBLE_MATRIX_3x3(mq2);
   TEST_DOUBLE_MATRIX_3x3(ml2);
   TEST_DOUBLE_MATRIX_3x3(mu2);
   TEST_DOUBLE_MATRIX_3x3(md2);
   TEST_DOUBLE_MATRIX_3x3(me2);
}

#define TEST_SINGLET(p)                                                 \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                            \
      BOOST_CHECK_EQUAL(mssm.get_##p(), 0.0);                           \
      Eigen::ArrayXd masses(mssm.get_DRbar_masses());                   \
      const double value = get_random();                                \
      masses(CMSSM_info::p) = value;                                    \
      mssm.set_DRbar_masses(masses);                                    \
      BOOST_CHECK_EQUAL(mssm.get_##p(), value);                         \
   }

#define TEST_MULTIPLET_2(p)                                             \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                            \
      for (int i = 0; i < 2; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), 0.0);                       \
      Eigen::ArrayXd masses(mssm.get_DRbar_masses());                   \
      Eigen::Array<double,2,1> p;                                       \
      for (int i = 0; i < 2; i++)                                       \
         p(i) = get_random();                                           \
      masses(CMSSM_info::p##_1) = p(0);                                 \
      masses(CMSSM_info::p##_2) = p(1);                                 \
      mssm.set_DRbar_masses(masses);                                    \
      for (int i = 0; i < 2; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), p(i));                      \
   }                                                                    \

#define TEST_MULTIPLET_3(p)                                             \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                            \
      for (int i = 0; i < 3; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), 0.0);                       \
      Eigen::ArrayXd masses(mssm.get_DRbar_masses());                   \
      Eigen::Array<double,3,1> p;                                       \
      for (int i = 0; i < 3; i++)                                       \
         p(i) = get_random();                                           \
      masses(CMSSM_info::p##_1) = p(0);                                 \
      masses(CMSSM_info::p##_2) = p(1);                                 \
      masses(CMSSM_info::p##_3) = p(2);                                 \
      mssm.set_DRbar_masses(masses);                                    \
      for (int i = 0; i < 3; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), p(i));                      \
   }                                                                    \

#define TEST_MULTIPLET_4(p)                                             \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                            \
      for (int i = 0; i < 4; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), 0.0);                       \
      Eigen::ArrayXd masses(mssm.get_DRbar_masses());                   \
      Eigen::Array<double,4,1> p;                                       \
      for (int i = 0; i < 4; i++)                                       \
         p(i) = get_random();                                           \
      masses(CMSSM_info::p##_1) = p(0);                                 \
      masses(CMSSM_info::p##_2) = p(1);                                 \
      masses(CMSSM_info::p##_3) = p(2);                                 \
      masses(CMSSM_info::p##_4) = p(3);                                 \
      mssm.set_DRbar_masses(masses);                                    \
      for (int i = 0; i < 4; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), p(i));                      \
   }                                                                    \

#define TEST_MULTIPLET_6(p)                                             \
   {                                                                    \
      CMSSM<Two_scale> mssm;                                            \
      for (int i = 0; i < 6; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), 0.0);                       \
      Eigen::ArrayXd masses(mssm.get_DRbar_masses());                   \
      Eigen::Array<double,6,1> p;                                       \
      for (int i = 0; i < 6; i++)                                       \
         p(i) = get_random();                                           \
      masses(CMSSM_info::p##_1) = p(0);                                 \
      masses(CMSSM_info::p##_2) = p(1);                                 \
      masses(CMSSM_info::p##_3) = p(2);                                 \
      masses(CMSSM_info::p##_4) = p(3);                                 \
      masses(CMSSM_info::p##_5) = p(4);                                 \
      masses(CMSSM_info::p##_6) = p(5);                                 \
      mssm.set_DRbar_masses(masses);                                    \
      for (int i = 0; i < 6; i++)                                       \
         BOOST_CHECK_EQUAL(mssm.get_##p(i), p(i));                      \
   }                                                                    \

BOOST_AUTO_TEST_CASE( test_CMSSM_mass_enum )
{
   // This test checks that the order of the enum entries
   // CMSSM_info::Masses matches the order of the particle masses in
   // the model class (as used in the set_DRbar_masses() and
   // get_DRbar_masses() functions).
   //
   // If this test works, then the user can fill a masses vector using
   // the enum entries and then set all model parameters at once via
   // the set() function.

   srand(1);

   {
      unsigned n_particles = 0;
      for (unsigned i = 0; i < CMSSM_info::NUMBER_OF_PARTICLES; i++)
         n_particles += CMSSM_info::particle_multiplicities[i];

      BOOST_CHECK_EQUAL(n_particles,
                        CMSSM<Two_scale>().get_DRbar_masses().rows());
      BOOST_CHECK_EQUAL(CMSSM_info::NUMBER_OF_MASSES,
                        CMSSM<Two_scale>().get_DRbar_masses().rows());
   }

   TEST_SINGLET(MVP);
   TEST_SINGLET(MVZ);
   TEST_SINGLET(MVWm);
   TEST_SINGLET(MVG);
   TEST_SINGLET(MGlu);
   TEST_MULTIPLET_2(Mhh);
   TEST_MULTIPLET_2(MAh);
   TEST_MULTIPLET_2(MHpm);
   TEST_MULTIPLET_3(MFu);
   TEST_MULTIPLET_3(MFd);
   TEST_MULTIPLET_3(MFe);
   TEST_MULTIPLET_3(MFv);
   TEST_MULTIPLET_4(MChi);
   TEST_MULTIPLET_2(MCha);
   TEST_MULTIPLET_6(MSu);
   TEST_MULTIPLET_6(MSd);
   TEST_MULTIPLET_6(MSe);
   TEST_MULTIPLET_3(MSv);
}
