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

#ifndef DEQ_SOLVER_H
#define DEQ_SOLVER_H

#include "linalg.h"

/**
 * @class Deq_solver
 * @author Alexander Voigt
 * @brief Differential equation solver interface
 */
class Deq_solver {
public:
   virtual ~Deq_solver() {}
   virtual DoubleVector solve(double) = 0;
};

#endif
