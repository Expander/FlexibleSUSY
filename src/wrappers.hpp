
#ifndef WRAPPERS_HPP
#define WRAPPERS_HPP

#include <complex>
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

inline double Abs(const std::complex<double>& z)
{
   return std::abs(z);
}

inline double AbsSqr(double z)
{
   return z * z;
}

inline double AbsSqr(const std::complex<double>& z)
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

inline std::complex<double> Conj(const std::complex<double>& a)
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

template <class Derived>
typename Eigen::MatrixBase<Derived>::PlainObject Diag(const Eigen::MatrixBase<Derived>& m)
{
   typename Eigen::MatrixBase<Derived>::PlainObject diag(m);
   for (int i = 0; i < m.rows(); ++i) {
      for (int k = 0; k < m.cols(); ++k) {
         if (i != k)
            diag(i,k) = 0.0;
      }
   }
   return diag;
}

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

inline double Re(const std::complex<double>& x)
{
   return std::real(x);
}

inline double Im(double x)
{
   return x;
}

inline double Im(const std::complex<double>& x)
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

#define UNITMATRIX(rows) Eigen::Matrix<double,rows,rows>::Identity()

inline double ZeroSqrt(double x)
{
   return (x > 0.0 ? std::sqrt(x) : 0.0);
}

namespace {
  inline double ZeroSqrt_d(double x)
  {
    return ZeroSqrt(x);
  }
}

template <typename Derived>
Derived ZeroSqrt(const Eigen::ArrayBase<Derived>& m)
{
   return m.unaryExpr(std::ptr_fun(ZeroSqrt_d));
}

}

#endif
