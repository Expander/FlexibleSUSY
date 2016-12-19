
/** \file numerics.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief numerical routines - differential equation solver, differentiator
   and function minimiser for instance
*/

#ifndef NUMERICS_H
#define NUMERICS_H

#include "utils.h"
#include "mycomplex.h"
#include <iostream>
#include <cmath>
#include "def.h"
#include "linalg.h"

/// func is user-supplied, h is an estimate of what step-size to start with
/// and err returns error flags
double calcDerivative(double (*func)(double), 
		     double x, double h, double *err);

/// f is user-defined function, minimum value returned in xmin. Based on a
/// golden section search
double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin);

/// Passarino-Veltman function definition
double b0(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b1(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b22(double p,  double m1, double m2, double q);
/// Passarino-Veltman function definition
double c0(double m1, double m2, double m3);
/// Passarino-Veltman function definition
double d27(double m1, double m2, double m3, double m4);
/// Passarino-Veltman function definition
double d0(double m1, double m2, double m3, double m4);

// inlined PV functions
inline double a0(double m, double q) {
   using std::fabs;
   using std::log;
  if (fabs(m) < softsusy::EPSTOL) return 0.;
  return sqr(m) * (1.0 - 2. * log(abs(m / q)));
}

inline double ffn(double p, double m1, double m2, double q) {
  return a0(m1, q) - 2.0 * a0(m2, q) - 
    (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) * 
    b0(p, m1, m2, q);
}

inline double gfn(double p, double m1, double m2, double q) {
  return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q) 
    - a0(m2, q); 
}

inline double hfn(double p, double m1, double m2, double q) {
  return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

inline double b22bar(double p, double m1, double m2, double q) {
  return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

double dilog(double x);

/// These three functions are for the calculation of 2-loop log pieces of g-2
/// of the muon
double fps(double z);
double fs(double z);
double ffbar(double z);

softsusy::Complex dilog(const softsusy::Complex & x);

/// Returns a 3 by 3 real mixing matrix. Input angles are standard CKM
/// parameterisation. If the phase d is not zero, the result is only an
/// approximation to the full complex matrix: see SOFTSUSY manual for
/// details. 
softsusy::DoubleMatrix display3x3RealMixing(double theta12, double theta13, 
				  double theta23, double d);

/// useful for 2-loop mb/mt corrections
double fin(double mm1, double mm2);
double den(double a, int b); /// 1/a^b
/// Given r and qt, carry out a Jacobi rotation on rows i and i+1 of each
/// matrix. a and b are parameters of the rotation: $\cos
/// \theta=a/\sqrt{a^2+b^2}, \sin \theta = b / \sqrt{a^2 + b^2}$.
void rotate(softsusy::DoubleMatrix & r, softsusy::DoubleMatrix & qt, int n, int i, float a, 
	    float b);

#endif

