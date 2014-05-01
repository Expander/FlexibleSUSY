#ifndef HIGGS_SCAN_OPTIONS_HPP
#define HIGGS_SCAN_OPTIONS_HPP

#define TANB_START    0.
#define TANB_STOP     50.
#define TANB_NPOINTS  100

#define M0_START    0.
#define M0_STOP     10000.
#define M0_NPOINTS  100

#include <cstdio>
#include <cstdlib>

struct Options {
   Options() { reset(); }
   Options(int argc, const char* argv[]) { reset(); parse(argc, argv); }
   ~Options() {}

   void parse(int, const char*[]);
   void reset();
   static bool starts_with(const std::string&, const std::string&);

   double lambda;
   double vs;
};

void Options::parse(int argc, const char* argv[])
{
   assert(argc > 0);

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);
      if (starts_with(option,"--lambda=")) {
         lambda = atof(option.substr(9).c_str());
      } else if (starts_with(option,"--vs=")) {
         vs = atof(option.substr(5).c_str());
      } else {
         ERROR("Unrecognized command line option: " << option);
         exit(EXIT_FAILURE);
      }
   }
}

void Options::reset()
{
   lambda = 0.1;
   vs = 10000.;
}

bool Options::starts_with(const std::string& str,
                          const std::string& prefix)
{
   return !str.compare(0, prefix.size(), prefix);
}

#endif
