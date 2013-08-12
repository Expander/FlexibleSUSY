
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_slha_output

#include <boost/test/unit_test.hpp>
#include <cstdlib>

int run_cmd(const std::string& cmd)
{
   if (!system(NULL)) {
      BOOST_MESSAGE("Error: command processor not available!");
      return -1;
   }

   const int status = system(cmd.c_str());

   BOOST_MESSAGE("Command \"" << cmd << "\" returned with exit code " << status);

   return status;
}

int run_point(const std::string& slha_file,
              const std::string& flexiblesusy_output_file,
              const std::string& softsusy_output_file)
{
   int status;

   status = run_cmd("./models/MSSM/run_MSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file=" + flexiblesusy_output_file);

   if (status) {
      BOOST_MESSAGE("FlexibleSUSY failed with exit code " << status);
      return status;
   }

   status = run_cmd("./examples/run_softpoint.x leshouches < " +
                    slha_file + " > " + softsusy_output_file);

   if (status) {
      BOOST_MESSAGE("Softsusy failed with exit code " << status);
      return status;
   }

   return 0;
}

BOOST_AUTO_TEST_CASE( test_slha_output )
{
   const std::string input_file("test/input_MSSM.slha2");
   const std::string output_file(input_file + ".spc");

   int status = run_point(input_file,
                          output_file + ".fs", output_file + ".ss");

   BOOST_REQUIRE_EQUAL(status, 0);
}
