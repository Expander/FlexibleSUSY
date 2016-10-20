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

#include "physical_input.hpp"

#include <cassert>

namespace flexiblesusy {

/**
 * Default constructor
 *
 * Calls reset() to initialize all physical input parameters to their
 * default values.
 */
Physical_input::Physical_input()
{
   reset();
}

double Physical_input::get(Input o) const
{
   return values.at(o);
}

Eigen::ArrayXd Physical_input::get() const
{
   Eigen::ArrayXd vec(values.size());

   std::copy(values.cbegin(), values.cend(), vec.data());

   return vec;
}

std::vector<std::string> Physical_input::get_names()
{
   std::vector<std::string> names(NUMBER_OF_INPUT_PARAMETERS);
   names[0] = "alpha_em(0)";
   names[1] = "mh_pole";
   return names;
}

void Physical_input::set(Input o, double value)
{
   values.at(o) = value;
}

void Physical_input::set(const Eigen::ArrayXd& vec)
{
   assert(vec.size() == values.size() && "Parameters array has wrong size");
   std::copy(vec.data(), vec.data() + vec.size(), values.begin());
}

/**
 * Resets all physical input parameters to their default values.
 *
 * | enum                             | possible values              | default value   |
 * |----------------------------------|------------------------------|-----------------|
 * | alpha_em_0                       | any positive double          | 1/137.035999074 |
 * | mh_pole                          | any positive double          | 125.09          |
 */
void Physical_input::reset()
{
   values[alpha_em_0] = 1./137.035999074;
   values[mh_pole] = 125.09;
}

} // namespace flexiblesusy
