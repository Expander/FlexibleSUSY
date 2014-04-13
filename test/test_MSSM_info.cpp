
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_info

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_model.hpp"
#include "MSSM_info.hpp"

using namespace flexiblesusy;

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

   // g1
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_g1(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::g1) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_g1(), value);
   }

   // g2
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_g2(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::g2) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_g2(), value);
   }

   // g3
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_g3(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::g3) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_g3(), value);
   }

   // Yu
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Yu(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> Yu;
      Yu << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::Yu00) = Yu(0,0);
      parameters(MSSM_info::Yu01) = Yu(0,1);
      parameters(MSSM_info::Yu02) = Yu(0,2);
      parameters(MSSM_info::Yu10) = Yu(1,0);
      parameters(MSSM_info::Yu11) = Yu(1,1);
      parameters(MSSM_info::Yu12) = Yu(1,2);
      parameters(MSSM_info::Yu20) = Yu(2,0);
      parameters(MSSM_info::Yu21) = Yu(2,1);
      parameters(MSSM_info::Yu22) = Yu(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Yu(i,k), Yu(i,k));
   }

   // Yd
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Yd(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> Yd;
      Yd << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::Yd00) = Yd(0,0);
      parameters(MSSM_info::Yd01) = Yd(0,1);
      parameters(MSSM_info::Yd02) = Yd(0,2);
      parameters(MSSM_info::Yd10) = Yd(1,0);
      parameters(MSSM_info::Yd11) = Yd(1,1);
      parameters(MSSM_info::Yd12) = Yd(1,2);
      parameters(MSSM_info::Yd20) = Yd(2,0);
      parameters(MSSM_info::Yd21) = Yd(2,1);
      parameters(MSSM_info::Yd22) = Yd(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Yd(i,k), Yd(i,k));
   }

   // Ye
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Ye(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> Ye;
      Ye << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::Ye00) = Ye(0,0);
      parameters(MSSM_info::Ye01) = Ye(0,1);
      parameters(MSSM_info::Ye02) = Ye(0,2);
      parameters(MSSM_info::Ye10) = Ye(1,0);
      parameters(MSSM_info::Ye11) = Ye(1,1);
      parameters(MSSM_info::Ye12) = Ye(1,2);
      parameters(MSSM_info::Ye20) = Ye(2,0);
      parameters(MSSM_info::Ye21) = Ye(2,1);
      parameters(MSSM_info::Ye22) = Ye(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_Ye(i,k), Ye(i,k));
   }

   // Mu
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_Mu(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::Mu) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_Mu(), value);
   }

   // vd
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_vd(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::vd) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_vd(), value);
   }

   // vu
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_vu(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::vu) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_vu(), value);
   }

   // ==================== soft parameters ====================

   // TYu
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYu(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> TYu;
      TYu << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::TYu00) = TYu(0,0);
      parameters(MSSM_info::TYu01) = TYu(0,1);
      parameters(MSSM_info::TYu02) = TYu(0,2);
      parameters(MSSM_info::TYu10) = TYu(1,0);
      parameters(MSSM_info::TYu11) = TYu(1,1);
      parameters(MSSM_info::TYu12) = TYu(1,2);
      parameters(MSSM_info::TYu20) = TYu(2,0);
      parameters(MSSM_info::TYu21) = TYu(2,1);
      parameters(MSSM_info::TYu22) = TYu(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYu(i,k), TYu(i,k));
   }

   // TYd
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYd(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> TYd;
      TYd << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::TYd00) = TYd(0,0);
      parameters(MSSM_info::TYd01) = TYd(0,1);
      parameters(MSSM_info::TYd02) = TYd(0,2);
      parameters(MSSM_info::TYd10) = TYd(1,0);
      parameters(MSSM_info::TYd11) = TYd(1,1);
      parameters(MSSM_info::TYd12) = TYd(1,2);
      parameters(MSSM_info::TYd20) = TYd(2,0);
      parameters(MSSM_info::TYd21) = TYd(2,1);
      parameters(MSSM_info::TYd22) = TYd(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYd(i,k), TYd(i,k));
   }

   // TYe
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYe(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> TYe;
      TYe << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::TYe00) = TYe(0,0);
      parameters(MSSM_info::TYe01) = TYe(0,1);
      parameters(MSSM_info::TYe02) = TYe(0,2);
      parameters(MSSM_info::TYe10) = TYe(1,0);
      parameters(MSSM_info::TYe11) = TYe(1,1);
      parameters(MSSM_info::TYe12) = TYe(1,2);
      parameters(MSSM_info::TYe20) = TYe(2,0);
      parameters(MSSM_info::TYe21) = TYe(2,1);
      parameters(MSSM_info::TYe22) = TYe(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_TYe(i,k), TYe(i,k));
   }

   // BMu
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_BMu(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::BMu) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_BMu(), value);
   }

   // MassB
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_MassB(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::MassB) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_MassB(), value);
   }

   // MassWB
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_MassWB(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::MassWB) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_MassWB(), value);
   }

   // MassG
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_MassG(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::MassG) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_MassG(), value);
   }

   // mHd2
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_mHd2(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::mHd2) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_mHd2(), value);
   }

   // mHu2
   {
      MSSM<Two_scale> mssm;
      BOOST_CHECK_EQUAL(mssm.get_mHu2(), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      const double value = 0.1;
      parameters(MSSM_info::mHu2) = value;
      mssm.set(parameters);

      BOOST_CHECK_EQUAL(mssm.get_mHu2(), value);
   }

   // mq2
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_mq2(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> mq2;
      mq2 << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::mq200) = mq2(0,0);
      parameters(MSSM_info::mq201) = mq2(0,1);
      parameters(MSSM_info::mq202) = mq2(0,2);
      parameters(MSSM_info::mq210) = mq2(1,0);
      parameters(MSSM_info::mq211) = mq2(1,1);
      parameters(MSSM_info::mq212) = mq2(1,2);
      parameters(MSSM_info::mq220) = mq2(2,0);
      parameters(MSSM_info::mq221) = mq2(2,1);
      parameters(MSSM_info::mq222) = mq2(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_mq2(i,k), mq2(i,k));
   }

   // ml2
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_ml2(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> ml2;
      ml2 << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::ml200) = ml2(0,0);
      parameters(MSSM_info::ml201) = ml2(0,1);
      parameters(MSSM_info::ml202) = ml2(0,2);
      parameters(MSSM_info::ml210) = ml2(1,0);
      parameters(MSSM_info::ml211) = ml2(1,1);
      parameters(MSSM_info::ml212) = ml2(1,2);
      parameters(MSSM_info::ml220) = ml2(2,0);
      parameters(MSSM_info::ml221) = ml2(2,1);
      parameters(MSSM_info::ml222) = ml2(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_ml2(i,k), ml2(i,k));
   }

   // mu2
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_mu2(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> mu2;
      mu2 << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::mu200) = mu2(0,0);
      parameters(MSSM_info::mu201) = mu2(0,1);
      parameters(MSSM_info::mu202) = mu2(0,2);
      parameters(MSSM_info::mu210) = mu2(1,0);
      parameters(MSSM_info::mu211) = mu2(1,1);
      parameters(MSSM_info::mu212) = mu2(1,2);
      parameters(MSSM_info::mu220) = mu2(2,0);
      parameters(MSSM_info::mu221) = mu2(2,1);
      parameters(MSSM_info::mu222) = mu2(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_mu2(i,k), mu2(i,k));
   }

   // md2
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_md2(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> md2;
      md2 << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::md200) = md2(0,0);
      parameters(MSSM_info::md201) = md2(0,1);
      parameters(MSSM_info::md202) = md2(0,2);
      parameters(MSSM_info::md210) = md2(1,0);
      parameters(MSSM_info::md211) = md2(1,1);
      parameters(MSSM_info::md212) = md2(1,2);
      parameters(MSSM_info::md220) = md2(2,0);
      parameters(MSSM_info::md221) = md2(2,1);
      parameters(MSSM_info::md222) = md2(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_md2(i,k), md2(i,k));
   }

   // me2
   {
      MSSM<Two_scale> mssm;
      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_me2(i,k), 0.0);

      Eigen::ArrayXd parameters(mssm.get_number_of_parameters());
      Eigen::Matrix<double,3,3> me2;
      me2 << 0.1, 0.2, 0.3,
         0.4, 0.5, 0.6,
         0.7, 0.8, 0.9;

      parameters(MSSM_info::me200) = me2(0,0);
      parameters(MSSM_info::me201) = me2(0,1);
      parameters(MSSM_info::me202) = me2(0,2);
      parameters(MSSM_info::me210) = me2(1,0);
      parameters(MSSM_info::me211) = me2(1,1);
      parameters(MSSM_info::me212) = me2(1,2);
      parameters(MSSM_info::me220) = me2(2,0);
      parameters(MSSM_info::me221) = me2(2,1);
      parameters(MSSM_info::me222) = me2(2,2);

      mssm.set(parameters);

      for (int i = 0; i < 3; i++)
         for (int k = 0; k < 3; k++)
            BOOST_CHECK_EQUAL(mssm.get_me2(i,k), me2(i,k));
   }
}
