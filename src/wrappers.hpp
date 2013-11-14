
#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "linalg.h"
#include <cmath>
#include <valarray>
#include <functional>
#include <Eigen/Core>

namespace flexiblesusy {

static const double Pi = M_PI;
static const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);

inline double Abs(double z)
{
   return std::fabs(z);
}

inline double Abs(const Complex& z)
{
   return std::abs(z);
}

inline double AbsSqr(double z)
{
   return z * z;
}

inline double AbsSqr(const Complex& z)
{
   return std::norm(z);
}

inline double AbsSqrt(double x)
{
   return std::sqrt(std::fabs(x));
}

inline double AbsSqrt_d(double x)
{
   return AbsSqrt(x);
}

inline DoubleVector AbsSqrt(const DoubleVector& x)
{
   return x.apply(AbsSqrt);
}

template <typename Derived>
Derived AbsSqrt(const Eigen::MatrixBase<Derived>& m)
{
   return m.unaryExpr(std::ptr_fun(AbsSqrt_d));
}

template <typename Derived>
Derived AbsSqrt(const Eigen::ArrayBase<Derived>& m)
{
   return m.unaryExpr(std::ptr_fun(AbsSqrt_d));
}

inline double ArcTan(double a)
{
   return atan(a);
}

inline double ArcSin(double a)
{
   return asin(a);
}

inline double ArcCos(double a)
{
   return acos(a);
}


inline double Conj(double a)
{
   return a;
}

inline Complex Conj(const Complex& a)
{
   return std::conj(a);
}

inline double Cos(double x)
{
   return cos(x);
}

inline double Sin(double x)
{
   return sin(x);
}

inline int Delta(int i, int j)
{
   return i == j;
}

inline int KroneckerDelta(int i, int j)
{
   return i == j;
}

Eigen::Matrix3d Diag(const Eigen::Matrix3d&);

inline double FiniteLog(double a)
{
   return a > std::numeric_limits<double>::epsilon() ? std::log(a) : 0;
}

inline double Log(double a)
{
   return std::log(a);
}

double MaxRelDiff(double, double);

template <class Derived>
double MaxRelDiff(const Eigen::MatrixBase<Derived>& a,
                  const Eigen::MatrixBase<Derived>& b)
{
   Derived sumTol;

   for (int i = 0; i < a.RowsAtCompileTime; i++) {
      const double max = maximum(a(i), b(i));
      if (std::fabs(max) > std::numeric_limits<double>::epsilon())
         sumTol(i) = fabs(1.0 - minimum(a(i), b(i)) / max);
      else
         sumTol(i) = 0.;
   }

   return sumTol.maxCoeff();
}

template <class Derived>
double MaxRelDiff(const Eigen::ArrayBase<Derived>& a,
                  const Eigen::ArrayBase<Derived>& b)
{
   Derived sumTol;

   for (int i = 0; i < a.RowsAtCompileTime; i++) {
      const double max = maximum(a(i), b(i));
      if (std::fabs(max) > std::numeric_limits<double>::epsilon())
         sumTol(i) = fabs(1.0 - minimum(a(i), b(i)) / max);
      else
         sumTol(i) = 0.;
   }

   return sumTol.maxCoeff();
}

template <typename Base, typename Exponent>
double Power(Base base, Exponent exp)
{
   return std::pow(base, exp);
}


inline double Re(double x)
{
   return x;
}

inline double Re(const Complex& x)
{
   return std::real(x);
}

inline double Im(double x)
{
   return x;
}

inline double Im(const Complex& x)
{
   return std::imag(x);
}

inline double Sqrt(double a)
{
   return std::sqrt(a);
}

template <typename T>
T Sqr(T a)
{
   return a * a;
}

template <typename Derived>
void Symmetrize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(m.RowsAtCompileTime == m.ColsAtCompileTime,
                 "Symmetrize is only defined for squared matrices");

   for (int i = 0; i < m.RowsAtCompileTime; i++)
      for (int k = 0; k < i; k++)
         m(i,k) = m(k,i);
}

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

inline ComplexMatrix Transpose(const ComplexMatrix& m)
{
   return m.transpose();
}

#define UNITMATRIX(rows) Eigen::Matrix<double,rows,rows>::Identity()

inline double ZeroSqrt(double x)
{
   return (x > 0.0 ? std::sqrt(x) : 0.0);
}

}

#endif
