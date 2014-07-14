
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_slha_io

#include <boost/test/unit_test.hpp>

#include "slha_io.hpp"
#include "linalg2.hpp"
#include "stopwatch.hpp"
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

/**
 * Creates a SLHA_io object with one block named `TestBlock' with
 *  `number_of_entries' entries of the form
 *
 *  Block TestBlock
 *     0  0
 *     1  1
 *     2  2
 *     3  3
 *
 * @param number_of_entries number of block entries
 *
 * @return a SLHA_io object with one block named `TestBlock'
 */
SLHA_io create_block(int number_of_entries)
{
   SLHAea::Coll coll;
   SLHAea::Block block;
   SLHA_io reader;
   std::string str = "Block TestBlock\n";

   for (int i = 0; i < number_of_entries; i++) {
      const std::string num(std::to_string(i));
      str += "   " + num + "  " + num + "\n";
   }

   block.str(str);
   coll.push_back(block);
   reader.set_data(coll);

   return reader;
}

void process_tuple(double* array, int key, double value)
{
   array[key] = value;
}

void check_array(double* array, int number_of_entries)
{
   for (int i = 0; i < number_of_entries; i++)
      BOOST_CHECK_EQUAL(array[i], i);
}

void reset_array(double* array, int number_of_entries, double value = 0.)
{
   for (int i = 0; i < number_of_entries; i++)
      array[i] = value;
}

/**
 * This test compares the speed of reading an array from a long SLHA
 * block using a) a tuple processor and b) a for loop.  The tuple
 * processor is more complicated to use, but is more than 1000 times
 * faster than the for loop.
 */
BOOST_AUTO_TEST_CASE( test_processor_vs_loop )
{
   using namespace std::placeholders;
   Stopwatch timer;
   double processor_time = 0., loop_time = 0.;
   const int number_of_entries = 10000;
   SLHA_io reader = create_block(number_of_entries);

   double array[number_of_entries];
   reset_array(array, number_of_entries);

   SLHA_io::Tuple_processor processor
      = std::bind(&process_tuple, array, _1, _2);

   // measure time to read the array using a tuple processor
   timer.start();
   reader.read_block("TestBlock", processor);
   timer.stop();
   processor_time = timer.get_time_in_seconds();

   check_array(array, number_of_entries);
   reset_array(array, number_of_entries);

   // measure time to read the array using a for loop
   timer.start();
   for (int i = 0; i < number_of_entries; i++) {
      array[i] = reader.read_entry("TestBlock", i);
   }
   timer.stop();
   loop_time = timer.get_time_in_seconds();

   check_array(array, number_of_entries);

   BOOST_MESSAGE("time using the tuple processor: " << processor_time << " s");
   BOOST_MESSAGE("time using the for loop: " << loop_time << " s");

   BOOST_CHECK_LT(1000. * processor_time, loop_time);
}

BOOST_AUTO_TEST_CASE( test_slha_mixing_matrix_convention )
{
   Eigen::Matrix<double, 2, 2> mass_matrix;
   Eigen::Matrix<std::complex<double>, 2, 2> Z;
   Eigen::Array<double, 2, 1> eigenvalues;

   mass_matrix(0,0) = 0;
   mass_matrix(0,1) = 1;
   mass_matrix(1,0) = 1;
   mass_matrix(1,1) = 0;

   // diagonalize with SARAH convention
   fs_diagonalize_symmetric(mass_matrix, eigenvalues, Z);

   BOOST_CHECK_GT(eigenvalues(0), 0.);
   BOOST_CHECK_GT(eigenvalues(1), 0.);

   BOOST_CHECK_GT(Z.imag().cwiseAbs().maxCoeff(), 0.);

   BOOST_CHECK((Z.row(0).imag().cwiseAbs().maxCoeff() > 0. &&
                Z.row(1).imag().cwiseAbs().maxCoeff() == 0.) ||
               (Z.row(1).imag().cwiseAbs().maxCoeff() > 0. &&
                Z.row(0).imag().cwiseAbs().maxCoeff() == 0.));

   // convert to SLHA convention
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(eigenvalues, Z);

   BOOST_CHECK(eigenvalues(0) < 0. || eigenvalues(1) < 0.);
   BOOST_CHECK_EQUAL(Z.imag().cwiseAbs().maxCoeff(), 0.);

   // reconstruct mass matrix
   Eigen::Matrix<std::complex<double>, 2, 2> reconstructed_mass_matrix;
   reconstructed_mass_matrix = Z.transpose() * eigenvalues.matrix().asDiagonal() * Z;

   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,0)), mass_matrix(0,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,1)), mass_matrix(0,1), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,0)), mass_matrix(1,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,1)), mass_matrix(1,1), 1.0e-15);

   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,1)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,1)), 0.);
}
