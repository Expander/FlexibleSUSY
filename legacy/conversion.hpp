
#ifndef CONVERSION_HPP
#define CONVERSION_HPP

#include "linalg.h"
#include <Eigen/Core>

namespace flexiblesusy {

Eigen::ArrayXd ToEigenArray(const DoubleVector&);
Eigen::ArrayXd ToEigenArray(double);
std::valarray<double> ToValarray(const DoubleVector&);
std::valarray<double> ToValarray(double);
Eigen::MatrixXd ToEigenMatrix(const DoubleMatrix&);

template<class Derived>
DoubleVector ToDoubleVector(const Eigen::ArrayBase<Derived>& a)
{
   DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

template<class Derived>
ComplexMatrix ToComplexMatrix(const Eigen::MatrixBase<Derived>& m)
{
   const int r = m.rows();
   const int c = m.cols();
   ComplexMatrix result(r,c);

   for (int i = 0; i < r; i++)
      for (int k = 0; k < c; k++)
         result(i+1, k+1) = m(i,k);

   return result;
}

template<class Derived>
DoubleMatrix ToDoubleMatrix(const Eigen::MatrixBase<Derived>& m)
{
   const int r = m.rows();
   const int c = m.cols();
   DoubleMatrix result(r,c);

   for (int i = 0; i < r; i++)
      for (int k = 0; k < c; k++)
         result(i+1, k+1) = m(i,k);

   return result;
}

} // namespace flexiblesusy

#endif
