
#include <iostream>
#include <fstream>
#include <string>

#include "test.h"
#include "slhaea.h"
#include "stopwatch.hpp"
#include "logger.hpp"

int run_cmd(const std::string& cmd)
{
   if (!system(NULL)) {
      ERROR("Error: command processor not available!");
      return -1;
   }

   const int status = system(cmd.c_str());

   if (status) {
      VERBOSE_MSG("Command \"" << cmd << "\" returned with exit code "
                  << status);
   }

   return status;
}

int run_point(const std::string& slha_file, double& fs_time, double& ss_time)
{
   int status;
   flexiblesusy::Stopwatch stopwatch;

   stopwatch.start();
   status = run_cmd("./models/MSSM/run_MSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file="
                    "test/test_benchmark.out.spc > /dev/null");
   stopwatch.stop();
   fs_time = stopwatch.get_time_in_seconds();

   if (status) {
      VERBOSE_MSG("FlexibleSUSY returned exit code " << status);
   }

   stopwatch.start();
   status = run_cmd("./examples/run_softpoint.x leshouches < " +
                    slha_file + " > test/test_benchmark.out.spc");
   stopwatch.stop();
   ss_time = stopwatch.get_time_in_seconds();

   if (status) {
      VERBOSE_MSG("Softsusy returned exit code " << status);
   }

   return 0;
}

SLHAea::Coll create_point(double tanBeta)
{
   std::ifstream ifs("test/test_benchmark.in.spc.in");
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

void test_tanbeta_scan()
{
   const unsigned num_points = 100;
   const double tanBeta_start = 2.;
   const double tanBeta_stop = 80.;
   const double tanBeta_step = (tanBeta_stop - tanBeta_start) / num_points;

   double fs_time_sum = 0., ss_time_sum = 0.;

   for (unsigned i = 0; i < num_points; i++) {
      const double tanBeta = tanBeta_start + i * tanBeta_step;
      const SLHAea::Coll coll(create_point(tanBeta));
      const std::string input_file("test/test_benchmark.in.spc");

      std::ofstream ofs(input_file);
      ofs << coll;
      ofs.close();

      double fs_time = 0., ss_time = 0.;

      INFO(">>> running point tan(beta) = " << tanBeta);
      run_point(input_file, fs_time, ss_time);

      INFO("\ttime: SS = " << ss_time << "s, FS = " << fs_time << "s"
           << ", rel. diff. = "
           << (fs_time - ss_time) * 100. / ss_time << "%");

      fs_time_sum += fs_time;
      ss_time_sum += ss_time;
   }

   fs_time_sum /= num_points;
   ss_time_sum /= num_points;

   INFO("Summary: average times (in seconds) \n"
        "  FlexibleSUSY: " << fs_time_sum << '\n' <<
        "  Softsusy    : " << ss_time_sum);

   TEST_GREATER(ss_time_sum, 2. * fs_time_sum);
}

int main()
{
   test_tanbeta_scan();

   return flexiblesusy::get_errors();
}
