
#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "linalg.h"
#include <cmath>
#include <valarray>
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

inline DoubleVector AbsSqrt(const DoubleVector& x)
{
   return x.apply(AbsSqrt);
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

Eigen::Matrix3d Diag(const Eigen::Matrix3d&);

void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);

// SVD
void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);

inline double FiniteLog(double a)
{
   return a > std::numeric_limits<double>::epsilon() ? std::log(a) : 0;
}

inline double Log(double a)
{
   return std::log(a);
}

double MaxRelDiff(double, double);
double MaxRelDiff(const DoubleVector&, const DoubleVector&);

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

inline DoubleMatrix Re(const DoubleMatrix& m)
{
   return m;
}

inline DoubleMatrix Re(const ComplexMatrix& m)
{
   return m.real();
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

void Symmetrize(DoubleMatrix&);

Eigen::ArrayXd ToEigenArray(const DoubleVector&);
Eigen::ArrayXd ToEigenArray(double);
std::valarray<double> ToValarray(const DoubleVector&);
std::valarray<double> ToValarray(double);
DoubleVector ToDoubleVector(const Eigen::ArrayXd&);
Eigen::MatrixXd ToEigenMatrix(const DoubleMatrix&);
DoubleMatrix ToDoubleMatrix(const Eigen::MatrixXd&);

inline DoubleMatrix Transpose(const DoubleMatrix& m)
{
   return m.transpose();
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

inline DoubleVector ZeroSqrt(const DoubleVector& x)
{
   return x.apply(ZeroSqrt);
}

}

#endif
