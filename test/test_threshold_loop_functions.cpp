// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_wrappers

#include <boost/test/unit_test.hpp>
#include "threshold_loop_functions.hpp"

using namespace flexiblesusy;
using namespace flexiblesusy::threshold_loop_functions;

BOOST_AUTO_TEST_CASE( test_C0 )
{
   struct Values {
      double m1, m2, m3;
   };

   Values value[] = {
      {0, 0, 0},
      {0, 0, 1},
      {0, 1, 1},
      {1, 1, 1},
      {1, 1, 0},
      {0, 1, 0},
      {1, 0, 1},
      {1, 2, 3}
   };

   for (int i = 0; i < sizeof(value)/sizeof(value[0]); i++) {
      BOOST_CHECK_EQUAL(C0(value[i].m1, value[i].m2, value[i].m3),
                        -Iabc(value[i].m1, value[i].m2, value[i].m3));
   }
}
