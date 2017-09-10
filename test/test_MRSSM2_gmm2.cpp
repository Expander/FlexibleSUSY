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
#define BOOST_TEST_MODULE test_MRSSM2_gmm2

#include <boost/test/unit_test.hpp>

#include "test_MRSSM2.hpp"

#include "wrappers.hpp"
#include "MRSSM2_a_muon.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu )
{
   MRSSM2_input_parameters input;

   // chargino dominance
   input.TanBeta = 10;
   input.Ms = 1000;
   input.LamTDInput = -1.0;
   input.LamTUInput = -1.0;
   input.LamSDInput = 1.1;
   input.LamSUInput = -1.1;
   input.MuDInput = 400;
   input.MuUInput = 400;
   input.BMuInput = Sqr(300);
   input.mq2Input << Sqr(1000), 0, 0,
		     0, Sqr(1000), 0,
		     0, 0, Sqr(1000);
   input.ml2Input << Sqr(500), 0, 0,
		     0, Sqr(500), 0,
		     0, 0, Sqr(500);
   input.md2Input << Sqr(1000), 0, 0,
		     0, Sqr(1000), 0,
		     0, 0, Sqr(1000);
   input.mu2Input << Sqr(1000), 0, 0,
		     0, Sqr(1000), 0,
		     0, 0, Sqr(1000);
   input.me2Input << Sqr(500), 0, 0,
		     0, Sqr(500), 0,
		     0, 0, Sqr(500);
   input.mS2Input = Sqr(2000);
   input.mT2Input = Sqr(3000);
   input.moc2Input = Sqr(1000);
   input.mRd2Input = Sqr(700);
   input.mRu2Input = Sqr(1000);
   input.MDBSInput = 1000;
   input.MDWBTInput = 500;
   input.MDGocInput = 1500;

   MRSSM2_slha<MRSSM2<Two_scale>> m = setup_MRSSM2(input);

   auto amu = MRSSM2_a_muon::calculate_a_muon(m);
   BOOST_CHECK_CLOSE(amu, -8.30e-11, 1.0);

   // neutralino dominance
   input.ml2Input << Sqr(8000), 0, 0,
		     0, Sqr(8000), 0,
		     0, 0, Sqr(8000);
   m = setup_MRSSM2(input);

   amu = MRSSM2_a_muon::calculate_a_muon(m);
   BOOST_CHECK_CLOSE(amu, 6.33e-12, 1.0);
}
