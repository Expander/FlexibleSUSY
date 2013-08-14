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

#include "program_options.hpp"

#include <cassert>

namespace flexiblesusy {

Program_options::Program_options()
   : values() // initializes all values to zero
{
   values[precision] = 1.0e-4;
   values[pole_mass_loop_order] = 1.;
}

double Program_options::get(Options o) const
{
   assert(o < NUMBER_OF_OPTIONS && "Option key out of range");
   return values[o];
}

void Program_options::set(Options o, double value)
{
   assert(o < NUMBER_OF_OPTIONS && "Option key out of range");
   values[o] = value;
}

} // namespace flexiblesusy
