
#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <cmath>
#include "linalg.h"
#include "numerics.h"
#include <Eigen/Core>

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

DoubleMatrix Diag(const DoubleMatrix&);
ComplexMatrix Diag(const ComplexMatrix&);
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
DoubleMatrix Re(const DoubleMatrix&);
DoubleMatrix Re(const ComplexMatrix&);

double Sqrt(double);

template <typename T>
T Sqr(T a)
{
   return a * a;
}

int ThetaStep(int, int);

DoubleMatrix Tp(const DoubleMatrix&);
ComplexMatrix Tp(const ComplexMatrix&);
DoubleMatrix Transpose(const DoubleMatrix&);
ComplexMatrix Transpose(const ComplexMatrix&);

double trace(const DoubleMatrix&);
Complex trace(const ComplexMatrix&);

#define UNITMATRIX(rows) Eigen::Matrix<double,rows,rows>::Identity()

double ZeroSqrt(double);

#endif
