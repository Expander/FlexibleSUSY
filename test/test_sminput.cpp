
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sminput

#include <boost/test/unit_test.hpp>

#include <fstream>
#include "slhaea.h"

BOOST_AUTO_TEST_CASE( test_reading_sminput )
{
   std::ifstream ifs("test/sminput.slha2");
   const SLHAea::Coll input(ifs);

   const double alpha_em = SLHAea::to<double>(input.at("SMINPUTS").at("1").at(1));
   const double G_mu     = SLHAea::to<double>(input.at("SMINPUTS").at("2").at(1));
   const double alpha_s  = SLHAea::to<double>(input.at("SMINPUTS").at("3").at(1));
   const double m_Z      = SLHAea::to<double>(input.at("SMINPUTS").at("4").at(1));
   const double m_b      = SLHAea::to<double>(input.at("SMINPUTS").at("5").at(1));
   const double m_t      = SLHAea::to<double>(input.at("SMINPUTS").at("6").at(1));
   const double m_tau    = SLHAea::to<double>(input.at("SMINPUTS").at("7").at(1));
   const double m_nu3    = SLHAea::to<double>(input.at("SMINPUTS").at("8").at(1));
   const double m_e      = SLHAea::to<double>(input.at("SMINPUTS").at("11").at(1));
   const double m_nu1    = SLHAea::to<double>(input.at("SMINPUTS").at("12").at(1));
   const double m_mu     = SLHAea::to<double>(input.at("SMINPUTS").at("13").at(1));
   const double m_nu2    = SLHAea::to<double>(input.at("SMINPUTS").at("14").at(1));
   const double m_d      = SLHAea::to<double>(input.at("SMINPUTS").at("21").at(1));
   const double m_u      = SLHAea::to<double>(input.at("SMINPUTS").at("22").at(1));
   const double m_s      = SLHAea::to<double>(input.at("SMINPUTS").at("23").at(1));
   const double m_c      = SLHAea::to<double>(input.at("SMINPUTS").at("24").at(1));

   BOOST_CHECK_EQUAL(alpha_em, 1.279340000e+02);
   BOOST_CHECK_EQUAL(G_mu    , 1.166370000e-05);
   BOOST_CHECK_EQUAL(alpha_s , 1.176000000e-01);
   BOOST_CHECK_EQUAL(m_Z     , 9.118760000e+01);
   BOOST_CHECK_EQUAL(m_b     , 4.200000000e+00);
   BOOST_CHECK_EQUAL(m_t     , 1.733000000e+02);
   BOOST_CHECK_EQUAL(m_tau   , 1.777000000e+00);
   BOOST_CHECK_EQUAL(m_nu3   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_e     , 5.109989020e-04);
   BOOST_CHECK_EQUAL(m_nu1   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_mu    , 1.056583570e-01);
   BOOST_CHECK_EQUAL(m_nu2   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_d     , 4.750000000e-03);
   BOOST_CHECK_EQUAL(m_u     , 2.400000000e-03);
   BOOST_CHECK_EQUAL(m_s     , 1.040000000e-01);
   BOOST_CHECK_EQUAL(m_c     , 1.270000000e+00);
}
