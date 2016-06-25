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

/********************* put types *********************/

/* put fundamental types */

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

#define MLPutRealMatrix(link,v,dim1,dim2)                  \
  do {                                                     \
    long dims[] = { dim1, dim2 };                          \
    MLPutDoubleArray(link, v, dims, NULL, NELEMS(dims));   \
  } while (0)

/* put real eigen types */

#define MLPutRealEigenArray(link,v,dim)                 \
   do {                                                 \
      double v_[dim];                                   \
      for (unsigned i = 0; i < dim; i++)                \
         v_[i] = v(i);                                  \
      MLPutRealList(link, v_, dim);                     \
   } while (0)

#define MLPutRealEigenVector(link,v,dim)                \
   MLPutRealEigenArray(link,v,dim)

#define MLPutRealEigenMatrix(link,M,dim1,dim2)               \
   do {                                                      \
      double M_[dim1][dim2];                                 \
      for (unsigned i = 0; i < dim1; i++)                    \
         for (unsigned k = 0; k < dim2; k++)                 \
            M_[i][k] = M(i, k);                              \
      MLPutRealMatrix(link, (double*)M_, dim1, dim2);        \
   } while (0)

/* put complex eigen types */

#define MLPutComplexEigenArray(link,v,dim)                       \
   do {                                                          \
      MLPutFunction(link, "List", dim);                          \
      for (unsigned i = 0; i < dim; i++)                         \
         MLPutComplex(link, std::real(v(i)), std::imag(v(i)));   \
   } while (0)

#define MLPutComplexEigenVector(link,v,dim)                \
   MLPutComplexEigenArray(link,v,dim)

#define MLPutComplexEigenMatrix(link,M,dim1,dim2)                       \
   do {                                                                 \
      MLPutFunction(link, "List", dim1);                                \
      for (unsigned i = 0; i < dim1; i++) {                             \
         MLPutFunction(link, "List", dim2);                             \
         for (unsigned k = 0; k < dim2; k++)                            \
            MLPutComplex(link, std::real(M(i,k)), std::imag(M(i,k)));   \
      }                                                                 \
   } while (0)

/********************* put rules *********************/

/* rules to fundamental types */

#define MLPutRule(link,name)                    \
   MLPutFunction(link, "Rule", 2);              \
   MLPutSymbol(link, (name))

#define MLPutRuleToReal(link,v,name) \
   MLPutRule(link, (name));          \
   MLPutReal(link, (v))

#define MLPutRuleToComplex(link,v,name) \
   MLPutRule(link, (name));             \
   MLPutComplex(link, std::real(v), std::imag(v))

#define MLPutRuleToInteger(link,v,name)         \
   MLPutRule(link, (name));                     \
   MLPutInteger(link, (v))

#define MLPutRuleToRealList(link,v,name,dim)    \
   MLPutRule(link, name); \
   MLPutRealList(link, v, dim)

#define MLPutRuleToRealMatrix(link,v,name,dim1, dim2)  \
   MLPutRule(link, name); \
   MLPutRealMatrix(link,v,dim1,dim2);

/* rules to Eigent types */

#define MLPutRuleToRealEigenArray(link,v,name,dim)      \
   MLPutRule(link, name);                               \
   MLPutRealEigenArray(link,v,dim)

#define MLPutRuleToRealEigenVector(link,v,name,dim)     \
   MLPutRule(link, name);                               \
   MLPutRealEigenVector(link,v,dim)

#define MLPutRuleToRealEigenMatrix(link,v,name,dim1,dim2)       \
   MLPutRule(link, name);                                       \
   MLPutRealEigenMatrix(link,v,dim1,dim2)

#define MLPutRuleToComplexEigenArray(link,v,name,dim)   \
   MLPutRule(link, name);                               \
   MLPutComplexEigenArray(link,v,dim)

#define MLPutRuleToComplexEigenVector(link,v,name,dim)  \
   MLPutRule(link, name);                               \
   MLPutComplexEigenVector(link,v,dim)

#define MLPutRuleToComplexEigenMatrix(link,v,name,dim1,dim2)    \
   MLPutRule(link, name);                                       \
   MLPutComplexEigenMatrix(link,v,dim1,dim2)

/*********** put rules to symbol with an additional head ***********/

/* rules to fundamental types */

#define MLPutRuleToHead(link,name,head)         \
   MLPutFunction(link, "Rule", 2);              \
   MLPutFunction(link, (head), 1);              \
   MLPutSymbol(link, (name))

#define MLPutRuleToHeadReal(link,v,name,head)   \
   MLPutRuleToHead(link, (name), (head));       \
   MLPutReal(link, (v))

#define MLPutRuleToHeadComplex(link,v,name,head)        \
   MLPutRuleToHead(link, (name), (head));               \
   MLPutComplex(link, std::real(v), std::imag(v))

#define MLPutRuleToHeadInteger(link,v,name,head)        \
   MLPutRuleToHead(link, (name), (head));               \
   MLPutInteger(link, (v))

#define MLPutRuleToHeadRealList(link,v,name,dim,head)   \
   MLPutRuleToHead(link, (name), (head));               \
   MLPutRealList(link, v, dim)

#define MLPutRuleToHeadRealMatrix(link,v,name,dim1,dim2,head)   \
   MLPutRuleToHead(link, (name), (head));                       \
   MLPutRealMatrix(link,v,dim1,dim2);

/* rules to Eigent types */

#define MLPutRuleToHeadRealEigenArray(link,v,name,dim,head)     \
   MLPutRuleToHead(link, (name), (head));                       \
   MLPutRealEigenArray(link,v,dim)

#define MLPutRuleToHeadRealEigenVector(link,v,name,dim,head)    \
   MLPutRuleToHead(link, (name), (head));                       \
   MLPutRealEigenVector(link,v,dim)

#define MLPutRuleToHeadRealEigenMatrix(link,v,name,dim1,dim2,head)      \
   MLPutRuleToHead(link, (name), (head));                               \
   MLPutRealEigenMatrix(link,v,dim1,dim2)

#define MLPutRuleToHeadComplexEigenArray(link,v,name,dim,head)  \
   MLPutRuleToHead(link, (name), (head));                       \
   MLPutComplexEigenArray(link,v,dim)

#define MLPutRuleToHeadComplexEigenVector(link,v,name,dim,head) \
   MLPutRuleToHead(link, (name), (head));                       \
   MLPutComplexEigenVector(link,v,dim)

#define MLPutRuleToHeadComplexEigenMatrix(link,v,name,dim1,dim2,head)   \
   MLPutRuleToHead(link, (name), (head));                               \
   MLPutComplexEigenMatrix(link,v,dim1,dim2)
