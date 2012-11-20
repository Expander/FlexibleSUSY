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

#ifndef TWO_SCALE_CONSTRAINT_H
#define TWO_SCALE_CONSTRAINT_H

#include "constraint.hpp"

class Two_scale;

template<>
class Constraint<Two_scale> {
public:
   virtual ~Constraint() {}
   virtual void apply() = 0;
   virtual double get_scale() const = 0;
   virtual void update_scale() = 0;
};

#endif
