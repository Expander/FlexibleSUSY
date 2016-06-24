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

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))

#define MLPutRule(link,name)                    \
   MLPutFunction(link, "Rule", 2);              \
   MLPutSymbol(link, (name))

#define MLPutRuleToReal(link,v,name) \
   MLPutRule(link, (name));          \
   MLPutReal(link, (v))

#define MLPutRuleToInteger(link,v,name)         \
   MLPutRule(link, (name));                     \
   MLPutInteger(link, (v))

#define MLPutRealMatrix(link,v,dim1,dim2)                  \
  do {                                                     \
    long dims[] = { dim1, dim2 };                          \
    MLPutDoubleArray(link, v, dims, NULL, NELEMS(dims));   \
  } while (0)

#define MLPutRuleToRealList(link,v,name,dim)    \
   MLPutFunction(link, "Rule", 2);              \
   MLPutSymbol(link, name);                     \
   MLPutRealList(link, v, dim)

#define MLPutRuleToRealMatrix(link,v,name,dim1, dim2)  \
   MLPutFunction(link, "Rule", 2);                     \
   MLPutSymbol(link, name);                            \
   MLPutRealMatrix(link,v,dim1,dim2);

#define MLPutRealEigenArray(link,v,dim)                 \
   do {                                                 \
      double v_[dim];                                   \
      for (unsigned i = 0; i < dim; i++)                \
         v_[i] = v(i);                                  \
      vLPutRealList(link, v_, dim);                     \
   } while (0)

#define MLPutRealEigenMatrix(link,M,dim1,dim2)               \
   do {                                                      \
      double M_[dim1][dim2];                                 \
      for (unsigned i = 0; i < dim1; i++)                    \
         for (unsigned k = 0; k < dim2; k++)                 \
            M_[i][k] = M(i, k);                              \
      MLPutRealMatrix(link, (double*)M_, dim1, dim2);        \
   } while (0)

#define MLPutComplexEigenMatrix(link,M,dim1,dim2)                       \
   do {                                                                 \
      MLPutFunction(link, "List", dim1);                                \
      for (unsigned i = 0; i < dim1; i++) {                             \
         MLPutFunction(link, "List", dim2);                             \
         for (unsigned k = 0; k < dim2; k++)                            \
            MLPutComplex(link, std::real(M(i,k)), std::imag(M(i,k)));   \
      }                                                                 \
   } while (0)

void MLPutComplex(MLINK link, double re, double im)
{
   if (im == 0.) {
      MLPutReal(link, re);
   } else {
      MLPutFunction(link, "Complex", 2);
      MLPutReal(link, re);
      MLPutReal(link, im);
   }
}
