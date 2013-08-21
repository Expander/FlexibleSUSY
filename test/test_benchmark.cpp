
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_benchmark

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>
#include <string>

#include "slhaea.h"
#include "stopwatch.hpp"

int run_cmd(const std::string& cmd)
{
   if (!system(NULL)) {
      BOOST_MESSAGE("Error: command processor not available!");
      return -1;
   }

   const int status = system(cmd.c_str());
   BOOST_REQUIRE(status == 0);

   if (status) {
      BOOST_MESSAGE("Command \"" << cmd << "\" returned with exit code "
                    << status);
   }

   return status;
}

int run_point(const std::string& slha_file)
{
   int status;
   flexiblesusy::Stopwatch stopwatch;
   double fs_time = 0., ss_time = 0.;

   stopwatch.start();
   status = run_cmd("./models/MSSM/run_MSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file=out.spc > /dev/null");
   stopwatch.stop();
   fs_time = stopwatch.get_time_in_seconds();

   if (status) {
      BOOST_MESSAGE("FlexibleSUSY failed with exit code " << status);
      return status;
   }

   stopwatch.start();
   status = run_cmd("./examples/run_softpoint.x leshouches < " +
                    slha_file + " > out.spc");
   stopwatch.stop();

   if (status) {
      BOOST_MESSAGE("Softsusy failed with exit code " << status);
      return status;
   }
   ss_time = stopwatch.get_time_in_seconds();

   BOOST_MESSAGE("\ttime: SS = " << ss_time << "s, FS = " << fs_time << "s"
                 << ", rel. diff. = "
                 << (fs_time - ss_time) * 100. / ss_time << "%");

   return 0;
}

SLHAea::Coll create_point(double tanBeta)
{
   std::ifstream ifs("test/slha_generic.slha2");
   SLHAea::Coll coll(ifs);
   SLHAea::Block minpar;

   const std::string str(
      "Block MINPAR\n"
      "   1   1.250000000e+02   # m0\n"
      "   2   5.000000000e+02   # m12\n"
      "   3   " + std::to_string(tanBeta) + "   # TanBeta\n"
      "   4   1.000000000e+00   # sign(Mu)\n"
      "   5   0.000000000e+00   # A0\n");

   minpar.str(str);
   coll.push_back(minpar);

   return coll;
}

BOOST_AUTO_TEST_CASE( test_tanbeta_scan )
{
   const unsigned num_points = 100;
   const double tanBeta_start = 2.;
   const double tanBeta_stop = 80.;
   const double tanBeta_step = (tanBeta_stop - tanBeta_start) / num_points;

   for (unsigned i = 0; i < num_points; i++) {
      const double tanBeta = tanBeta_start + i * tanBeta_step;
      const SLHAea::Coll coll(create_point(tanBeta));
      const std::string input_file("test/tmp_point.slha2");

      std::ofstream ofs(input_file);
      ofs << coll;
      ofs.close();

      BOOST_MESSAGE(">>> running point tan(beta) = " << tanBeta);
      run_point(input_file);
   }
}
