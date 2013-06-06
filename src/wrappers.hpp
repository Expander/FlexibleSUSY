
#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <cmath>
#include "linalg.h"
#include "numerics.h"

static const double Pi = M_PI;

double AbsSqr(const Complex&);

double ArcTan(double);
double ArcSin(double);
double ArcCos(double);

DoubleMatrix Adj(const DoubleMatrix&);
ComplexMatrix Adj(const ComplexMatrix&);

double Conj(double);
Complex Conj(const Complex&);
DoubleMatrix Conj(const DoubleMatrix&);
ComplexMatrix Conj(const ComplexMatrix&);

double Cos(double);
double Sin(double);

int Delta(int, int);

void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);
void DiagonalizeUnsorted(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void DiagonalizeUnsorted(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);

// SVD
void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);

double Log(double);
double Mass2(double);

double MaxRelDiff(double, double);
double MaxRelDiff(const DoubleVector&, const DoubleVector&);

template <typename Base, typename Exponent>
double Power(Base base, Exponent exp)
{
   return std::pow(base, exp);
}

double Re(double);
double Re(const Complex&);

double Sqrt(double);

template <typename T>
T Sqr(T a)
{
   return a * a;
}

int ThetaStep(int, int);

DoubleMatrix Tp(const DoubleMatrix&);
ComplexMatrix Tp(const ComplexMatrix&);

double trace(const DoubleMatrix&);
Complex trace(const ComplexMatrix&);

#define UNITMATRIX(rows) \
   unitMatrix<rows>()

template <int rows>
DoubleMatrix unitMatrix()
{
   DoubleMatrix u(rows,rows);
   for (int i = 1; i <= rows; ++i)
      u(i,i) = 1.0;
   return u;
}

#endif
