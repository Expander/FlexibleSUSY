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

#include "gsl_vector.hpp"
#include <cassert>

namespace flexiblesusy {

GSL_vector::GSL_vector()
   : vec(NULL)
{
}

GSL_vector::GSL_vector(std::size_t size)
{
   if (size)
      vec = gsl_vector_calloc(size);
   else
      vec = NULL;
}

GSL_vector::GSL_vector(const GSL_vector& other)
{
   vec = gsl_vector_alloc(other.size());
   gsl_vector_memcpy(vec, other.vec);
}

GSL_vector::GSL_vector(GSL_vector&& other)
{
   vec = other.vec;
   other.vec = NULL;
}

GSL_vector::~GSL_vector()
{
   gsl_vector_free(vec);
}

const GSL_vector& GSL_vector::operator=(const GSL_vector& other)
{
   if (this != &other) {
      gsl_vector_free(vec);
      vec = gsl_vector_alloc(other.size());
      gsl_vector_memcpy(vec, other.vec);
   }

   return *this;
}

double& GSL_vector::operator[](std::size_t n)
{
   assert(vec);
   assert(n < vec->size && "GSL_vector::operator[]: index out of range");
   return *gsl_vector_ptr(vec, n);
}

double GSL_vector::operator[](std::size_t n) const
{
   assert(vec);
   assert(n < vec->size && "GSL_vector::operator[]: index out of range");
   return gsl_vector_get(vec, n);
}

std::size_t GSL_vector::size() const
{
   if (!vec) return 0;
   return vec->size;
}

const gsl_vector* GSL_vector::raw() const
{
   return vec;
}

gsl_vector* GSL_vector::raw()
{
   return vec;
}

std::ostream& operator<<(std::ostream& ostr, const GSL_vector& vec)
{
   std::cout << "(";

   for (std::size_t i = 0; i < vec.size(); i++) {
      std::cout << vec[i];
      if (i < vec.size() - 1)
         std::cout << ", ";
   }

   std::cout << ")";
}

} // namespace flexiblesusy
