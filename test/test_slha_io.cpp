
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_slha_io

#include <boost/test/unit_test.hpp>

#include "slha_io.hpp"
#include <Eigen/Core>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_read_entry )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block A\n"
      "   1   1.2      # comment 1\n"
      "   2   2.3E+2   # comment 2";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   BOOST_CHECK_EQUAL(reader.read_entry("A",1), 1.2);
   BOOST_CHECK_EQUAL(reader.read_entry("A",2), 230.);
}

BOOST_AUTO_TEST_CASE( test_read_block )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block Matrix\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::MatrixXd matrix(Eigen::MatrixXd::Zero(2,2));
   reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(matrix(0,0), 1.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 2.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 4.0);
}

BOOST_AUTO_TEST_CASE( test_read_block_dense )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block Matrix\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::Matrix<double,2,2> matrix(Eigen::Matrix<double,2,2>::Zero());
   reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(matrix(0,0), 1.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 2.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 4.0);
}
