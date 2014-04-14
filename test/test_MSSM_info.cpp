
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_info

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_model.hpp"
#include "MSSM_info.hpp"

using namespace flexiblesusy;

#define TEST_SCALAR(p)                                                  \
   {                                                                    \
      MSSM<Two_scale> mssm;                                             \
      BOOST_CHECK_EQUAL(mssm.get_##p(), 0.0);                           \
      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());       \
      const double value = 0.1;                                         \
      parameters(MSSM_info::p) = value;                                 \
      mssm.set(parameters);                                             \
      BOOST_CHECK_EQUAL(mssm.get_##p(), value);                         \
   }

#define TEST_DOUBLE_MATRIX_3x3(p)                                       \
   {                                                                    \
      MSSM<Two_scale> mssm;                                             \
      for (int i = 0; i < 3; i++)                                       \
         for (int k = 0; k < 3; k++)                                    \
            BOOST_CHECK_EQUAL(mssm.get_##p(i,k), 0.0);                  \
      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());       \
      Eigen::Matrix<double,3,3> p;                                      \
      p << 0.1, 0.2, 0.3,                                               \
         0.4, 0.5, 0.6,                                                 \
         0.7, 0.8, 0.9;                                                 \
      parameters(MSSM_info::p##00) = p(0,0);                            \
      parameters(MSSM_info::p##01) = p(0,1);                            \
      parameters(MSSM_info::p##02) = p(0,2);                            \
      parameters(MSSM_info::p##10) = p(1,0);                            \
      parameters(MSSM_info::p##11) = p(1,1);                            \
      parameters(MSSM_info::p##12) = p(1,2);                            \
      parameters(MSSM_info::p##20) = p(2,0);                            \
      parameters(MSSM_info::p##21) = p(2,1);                            \
      parameters(MSSM_info::p##22) = p(2,2);                            \
      mssm.set(parameters);                                             \
      for (int i = 0; i < 3; i++)                                       \
         for (int k = 0; k < 3; k++)                                    \
            BOOST_CHECK_EQUAL(mssm.get_##p(i,k), p(i,k));               \
   }                                                                    \

BOOST_AUTO_TEST_CASE( test_MSSM_parameter_enum )
{
   // This test checks that the order of the enum entries
   // MSSM_info::Parameters matches the order of the parameters in the
   // model class (as used in the set() and display() functions).
   //
   // If this test works, then the user can fill a parameter vector
   // using the enum entries and then set all model parameters at once
   // via the set() function.

   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(MSSM_info::NUMBER_OF_PARAMETERS,
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
