
#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "linalg.h"
#include <cmath>
#include <Eigen/Core>

static const double Pi = M_PI;

double AbsSqr(const Complex&);

double ArcTan(double);
double ArcSin(double);
double ArcCos(double);

double Conj(double);
Complex Conj(const Complex&);

double Cos(double);
double Sin(double);

int Delta(int, int);

Eigen::Matrix3d Diag(const Eigen::Matrix3d&);

void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, DoubleMatrix& , DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, DoubleVector&);

// SVD
void Diagonalize(const DoubleMatrix&, DoubleMatrix& , DoubleMatrix& , DoubleVector&);
void Diagonalize(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);
void Diagonalize2by2(const DoubleMatrix&, ComplexMatrix&, ComplexMatrix&, DoubleVector&);

double Log(double);

double MaxRelDiff(double, double);
double MaxRelDiff(const DoubleVector&, const DoubleVector&);

template <typename Base, typename Exponent>
double Power(Base base, Exponent exp)
{
   return std::pow(base, exp);
}

double Re(double);
double Re(const Complex&);
DoubleMatrix Re(const DoubleMatrix&);
DoubleMatrix Re(const ComplexMatrix&);

double Sqrt(double);

template <typename T>
T Sqr(T a)
{
   return a * a;
}

int ThetaStep(int, int);

DoubleMatrix Transpose(const DoubleMatrix&);
ComplexMatrix Transpose(const ComplexMatrix&);

#define UNITMATRIX(rows) Eigen::Matrix<double,rows,rows>::Identity()

double ZeroSqrt(double);

#endif
