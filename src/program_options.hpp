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
   /// Spectrum generator options
   enum Options : unsigned {
      precision,            ///< overall precision goal
      max_iterations,       ///< maximum number of iterations (0 = automatic)
      algorithm,            ///< RG solver algorithm (0 = two-scale)
      calculate_sm_masses,  ///< calculate Standard Model pole masses
      pole_mass_loop_order, ///< loop-order for calculation of pole masses
      ewsb_loop_order,      ///< loop-order for solving the EWSB eqs.
      NUMBER_OF_OPTIONS     ///< number of possible options
   };

   Program_options();
   ~Program_options() {}

   double get(Options) const;
   void set(Options, double);
   void reset();

private:
   double values[NUMBER_OF_OPTIONS];
};

} // namespace flexiblesusy

#endif
