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

#ifndef GSL_VECTOR_H
#define GSL_VECTOR_H

#include <gsl/gsl_vector.h>
#include <iostream>
#include <cstddef>

namespace flexiblesusy {

class GSL_vector {
public:
   GSL_vector();
   explicit GSL_vector(std::size_t);
   GSL_vector(const GSL_vector&);
   GSL_vector(GSL_vector&&) noexcept;
   ~GSL_vector();

   const GSL_vector& operator=(const GSL_vector&);
   GSL_vector& operator=(GSL_vector&&);
   double& operator[](std::size_t);      ///< element read/write access
   double operator[](std::size_t) const; ///< element read access

   void assign(const gsl_vector*); ///< assign from gsl_vector
   const gsl_vector* raw() const;  ///< get raw pointer
   gsl_vector* raw();              ///< get raw pointer
   void set_all(double);           ///< set all elemets to same value
   std::size_t size() const;

private:
   gsl_vector* vec;                ///< raw gsl_vector

   void move_assign(GSL_vector&&); ///< move assign
   void range_check(std::size_t) const;
};

std::ostream& operator<<(std::ostream&, const GSL_vector&);

} // namespace flexiblesusy

#endif
