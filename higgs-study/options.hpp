#ifndef HIGGS_SCAN_OPTIONS_HPP
#define HIGGS_SCAN_OPTIONS_HPP

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
   double kappa;
   double vs;

   double tanb_start, tanb_stop;
   double m0_start, m0_stop;
   unsigned tanb_npoints;
   unsigned m0_npoints;
};

void Options::parse(int argc, const char* argv[])
{
   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);
      if (starts_with(option,"--lambda=")) {
         lambda = atof(option.substr(9).c_str());
      } else if (starts_with(option,"--vs=")) {
         vs = atof(option.substr(5).c_str());
      } else if (starts_with(option,"--kappa=")) {
         kappa = atof(option.substr(8).c_str());
      } else if (starts_with(option,"--tanb-start=")) {
         tanb_start = atof(option.substr(13).c_str());
      } else if (starts_with(option,"--tanb-stop=")) {
         tanb_stop = atof(option.substr(12).c_str());
      } else if (starts_with(option,"--tanb-npoints=")) {
         tanb_npoints = atof(option.substr(15).c_str());
      } else if (starts_with(option,"--m0-start=")) {
         m0_start = atof(option.substr(11).c_str());
      } else if (starts_with(option,"--m0-stop=")) {
         m0_stop = atof(option.substr(10).c_str());
      } else if (starts_with(option,"--m0-npoints=")) {
         m0_npoints = atof(option.substr(13).c_str());
      } else {
         ERROR("Unrecognized command line option: " << option);
         exit(EXIT_FAILURE);
      }
   }
}

void Options::reset()
{
   lambda = 0.1;
   kappa = 0.1;
   vs = 10000.;

   tanb_start   = 0.;
   tanb_stop    = 50.;
   tanb_npoints = 10;

   m0_start     = 0.;
   m0_stop      = 10000.;
   m0_npoints   = 10;
}

bool Options::starts_with(const std::string& str,
                          const std::string& prefix)
{
   return !str.compare(0, prefix.size(), prefix);
}

#endif
