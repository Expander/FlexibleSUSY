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

/**
 * @file observables.hpp
 * @brief contains class Observables
 */

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

namespace flexiblesusy {

struct Observables {
   Observables();
   void clear(); ///< sets all observables to zero

   double amu;   ///< a_mu = (g-2)/2 of the muon
};

} // namespace flexiblesusy

#endif
