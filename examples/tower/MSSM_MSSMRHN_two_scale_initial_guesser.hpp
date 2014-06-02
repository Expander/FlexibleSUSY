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

#ifndef MSSM_MSSMRHN_TWO_SCALE_INITIAL_GUESSER_H
#define MSSM_MSSMRHN_TWO_SCALE_INITIAL_GUESSER_H

#include "MSSM_MSSMRHN_initial_guesser.hpp"
#include "MSSMRHN_input_parameters.hpp"
#include "MSSM_two_scale_low_scale_constraint.hpp"
#include "MSSM_two_scale_susy_scale_constraint.hpp"
#include "MSSMRHN_two_scale_high_scale_constraint.hpp"
#include "MSSM_MSSMRHN_two_scale_matching.hpp"
#include "two_scale_initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class MSSM;

template <class T>
class MSSMRHN;

class Two_scale;

template<>
class MSSM_MSSMRHN_initial_guesser<Two_scale> : public Initial_guesser<Two_scale> {
public:
   MSSM_MSSMRHN_initial_guesser(MSSM<Two_scale>*, MSSMRHN<Two_scale>*,
				const MSSMRHN_input_parameters&,
				const QedQcd&,
				const MSSM_low_scale_constraint<Two_scale>&,
				const MSSM_susy_scale_constraint<Two_scale>&,
				const MSSMRHN_high_scale_constraint<Two_scale>&,
				const MSSM_MSSMRHN_matching<Two_scale>&);
   virtual ~MSSM_MSSMRHN_initial_guesser();
   virtual void guess();

private:
   MSSM<Two_scale>* model_1;
   MSSMRHN<Two_scale>* model_2;
   MSSMRHN_input_parameters input_pars;
   QedQcd oneset;
   MSSM_low_scale_constraint<Two_scale> low_constraint_1;
   MSSM_susy_scale_constraint<Two_scale> susy_constraint_1;
   MSSMRHN_high_scale_constraint<Two_scale> high_constraint_2;
   MSSM_MSSMRHN_matching<Two_scale> matching;

   void guess_susy_parameters();
   void guess_soft_parameters();
};

} // namespace flexiblesusy

#endif
