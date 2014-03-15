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

/**
 * Default constructor
 *
 * Calls reset() to initialize all program options to their default
 * values.
 */
Program_options::Program_options()
{
   reset();
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

/**
 * Resets all program options to their defaults
 *
 * | enum                 | possible values              | default value   |
 * |----------------------|------------------------------|-----------------|
 * | precision            | any positive double          | 1.0e-4          |
 * | max_iterations       | any positive double          | 0 (= automatic) |
 * | algorithm            | 0 (two-scale) or 1 (lattice) | 0 (= two-scale) |
 * | calculate_sm_masses  | 0 or 1                       | 0 (= no)        |
 * | pole_mass_loop_order | 0, 1, 2                      | 2 (= 2-loop)    |
 * | ewsb_loop_order      | 0, 1, 2                      | 2 (= 2-loop)    |
 */
void Program_options::reset()
{
   values[precision]            = 1.0e-4;
   values[max_iterations]       = 0.; // 0 = automatic
   values[algorithm]            = 0.; // 0 = two-scale
   values[calculate_sm_masses]  = 0.; // 0 = false
   values[pole_mass_loop_order] = 2.;
   values[ewsb_loop_order]      = 2.;
}

} // namespace flexiblesusy
