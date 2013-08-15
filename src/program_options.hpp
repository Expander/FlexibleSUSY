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

#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H

namespace flexiblesusy {

/**
 * @class Program_options
 * @brief stores the program options
 *
 * This class stores all program options which can be changed via the
 * SLHA input file.
 */
class Program_options {
public:
   enum Options : unsigned { precision, max_iterations,
         algorithm, calculate_sm_masses, pole_mass_loop_order,
         ewsb_loop_order, NUMBER_OF_OPTIONS };

   Program_options();
   ~Program_options() {}

   double get(Options) const;
   void set(Options, double);

private:
   double values[NUMBER_OF_OPTIONS];
};

} // namespace flexiblesusy

#endif
