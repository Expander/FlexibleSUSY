
#include "MSSM_two_scale_susy_parameters.hpp"
#include "ew_input.hpp"
#include "logger.hpp"

using namespace flexiblesusy;

int main()
{
   MSSM_susy_parameters susy;
   susy.set_scale(Electroweak_constants::MZ);

   susy.run_to(1.0e16);
   susy.set_g1(0.1);
   susy.set_g2(0.1);
   susy.run_to(Electroweak_constants::MZ);

   INFO("g1(MZ) = " << susy.get_g1() << ", g2(MZ) = " << susy.get_g2());

   return 0;
}
