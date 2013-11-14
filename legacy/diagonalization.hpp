
#ifndef DIAGONALIZATION_HPP
#define DIAGONALIZATION_HPP

#include "linalg.h"

namespace flexiblesusy {

void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);

// SVD
void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);

double MaxRelDiff(const DoubleVector&, const DoubleVector&);

void Symmetrize(DoubleMatrix&);

inline DoubleMatrix Re(const DoubleMatrix& m)
{
   return m;
}

inline DoubleMatrix Re(const ComplexMatrix& m)
{
   return m.real();
}

inline DoubleMatrix Transpose(const DoubleMatrix& m)
{
   return m.transpose();
}

DoubleVector ZeroSqrt(const DoubleVector&);

} // namespace flexiblesusy

#endif
