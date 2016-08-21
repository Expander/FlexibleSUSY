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

namespace flexiblesusy {

class GSL_vector {
public:
   GSL_vector();
   explicit GSL_vector(std::size_t);
   GSL_vector(const GSL_vector&);
   GSL_vector(GSL_vector&&);
   ~GSL_vector();

   const GSL_vector& operator=(const GSL_vector&);
   double& operator[](std::size_t);
   double operator[](std::size_t) const;

   void assign(const gsl_vector*);
   bool all_finite() const;
   const gsl_vector* raw() const;
   gsl_vector* raw();
   void set_all(double);
   std::size_t size() const;

private:
   gsl_vector* vec;
};

std::ostream& operator<<(std::ostream&, const GSL_vector&);

} // namespace flexiblesusy

#endif
