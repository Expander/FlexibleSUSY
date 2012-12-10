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

#ifndef MSSM_TWO_SCALE_INITIAL_GUESSER_H
#define MSSM_TWO_SCALE_INITIAL_GUESSER_H

#include "two_scale_initial_guesser.hpp"
#include "lowe.h"
#include "linalg.h"

template<class T> class Mssm;
class Two_scale;

class Mssm_initial_guesser : public Initial_guesser<Two_scale> {
public:
   Mssm_initial_guesser(Mssm<Two_scale>*, const QedQcd&, double, double, int, const DoubleVector&);
   virtual ~Mssm_initial_guesser();
   virtual void guess();

private:
   Mssm<Two_scale>* mssm;     ///< Mssm model
   const QedQcd oneset;       ///< low-energy parameters
   double mxGuess;            ///< guessed GUT scale
   double tanb;               ///< tan(beta)
   int sgnMu;                 ///< sign of mu
   const DoubleVector pars;   ///< GUT parameters m0, m12, a0
};

#endif
