
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include "test.h"
#include "run_cmd.hpp"
#include "slhaea.h"
#include "stopwatch.hpp"
#include "logger.hpp"

struct Data {
   Data() : number_of_valid_points(0), sum_of_times(0.) {}

   double time;
   int error;
   unsigned number_of_valid_points;
   double sum_of_times;
};

bool is_valid_spectrum(const std::string& slha_file)
{
   std::ifstream ifs(slha_file.c_str());
   SLHAea::Coll coll(ifs);

   // find SPINFO block
   SLHAea::Coll::const_iterator block =
      coll.find(coll.cbegin(), coll.cend(), "SPINFO");

   if (block == coll.cend())
      throw std::string("Error: SPINFO block not found in file ") + slha_file;

   for (SLHAea::Block::const_iterator line = block->cbegin(),
           end = block->cend(); line != end; ++line) {
      if (line->is_data_line() && line->size() >= 2 &&
          (*line)[0] == "4" && (*line)[1] != "")
         return false;
   }

   return true;
}

int run_point(const std::string& slha_file, Data& fs_data, Data& ss_data)
{
   int status;
   flexiblesusy::Stopwatch stopwatch;

   const std::string slha_output_file("test/test_benchmark.out.spc");

   stopwatch.start();
   status = run_cmd("./models/CMSSM/run_CMSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file=" + slha_output_file +
                    " > /dev/null 2>&1");
   stopwatch.stop();

   fs_data.time = stopwatch.get_time_in_seconds();
   fs_data.error = status;

   if (!fs_data.error) {
      // look for errors in the SLHA output file
      fs_data.error = !is_valid_spectrum(slha_output_file);
      if (!fs_data.error) {
         fs_data.number_of_valid_points++;
         fs_data.sum_of_times += fs_data.time;
      }
   }

   stopwatch.start();
   status = run_cmd("./models/SoftsusyNMSSM/run_softpoint.x leshouches < " +
                    slha_file + " > " + slha_output_file);
   stopwatch.stop();

   ss_data.time = stopwatch.get_time_in_seconds();
   ss_data.error = status;

   if (!ss_data.error) {
      // look for errors in the SLHA output file
      ss_data.error = !is_valid_spectrum(slha_output_file);
      if (!ss_data.error) {
         ss_data.number_of_valid_points++;
         ss_data.sum_of_times += ss_data.time;
      }
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
      "   3   " + boost::lexical_cast<std::string>(tanBeta) + "   # TanBeta\n"
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

   Data fs_data, ss_data;

   printf("%10s %30s %30s \n", "tan(beta)",
          "Softsusy / s (status)",
          "FlexibleSUSY / s (status)");

   for (unsigned i = 0; i < num_points; i++) {
      const double tanBeta = tanBeta_start + i * tanBeta_step;
      const SLHAea::Coll coll(create_point(tanBeta));
      const std::string input_file("test/test_benchmark.in.spc");

      std::ofstream ofs(input_file);
      ofs << coll;
      ofs.close();

      run_point(input_file, fs_data, ss_data);

      printf("%10g %24g (%3d) %24g (%3d)\n", tanBeta,
             ss_data.time, ss_data.error,
             fs_data.time, fs_data.error);
   }

   const double fs_average_time = fs_data.sum_of_times / fs_data.number_of_valid_points;
   const double ss_average_time = ss_data.sum_of_times / ss_data.number_of_valid_points;

   INFO("Summary: average times (in seconds) \n"
        "  FlexibleSUSY: " << fs_average_time <<
        " (" << fs_data.number_of_valid_points << "/" << num_points << " points)\n" <<
        "  Softsusy    : " << ss_average_time <<
        " (" << ss_data.number_of_valid_points << "/" << num_points << " points)");

   TEST_GREATER(ss_average_time, 1.6 * fs_average_time);
}

int main()
{
   test_tanbeta_scan();

   return flexiblesusy::get_errors();
}
