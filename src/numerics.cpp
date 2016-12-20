
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

/// Comment if you want default softsusy behaviour
// #define USE_LOOPTOOLS

#include <Eigen/Dense>
#include "numerics.h"
#include "dilog.hpp"
#include "error.hpp"
#include "rk.hpp"
#ifdef USE_LOOPTOOLS
#include "clooptools.h"
#endif
#include <cmath>

using namespace softsusy;
using namespace Eigen;

namespace flexiblesusy {

class SoftsusyNumericsError : public Error {
public:
   explicit SoftsusyNumericsError(std::string msg_)
      : msg(msg_)
      {}
   virtual ~SoftsusyNumericsError() {}
   virtual std::string what() const {
      return msg;
   }
private:
   std::string msg;
};

} // namespace flexiblesusy

namespace {

// Integration routine needs these variables
static double mtInt, pInt, m1Int, m2Int;
static int nInt;

inline void shft2(double & a, double & b, double c) { a = b; b = c; }

inline void shft3(double & a, double & b, double & c, double d) { 
  a = b; b = c; c = d;
}

double refnfn(double x, double p, double m1, double m2, double q) noexcept
{
  const static std::complex<double> iEpsilon(0.0, TOLERANCE * 1.0e-20);

  return std::real(x *
    std::log(((1 - x) * sqr(m1) + x * sqr(m2)
              - x * (1 - x) * sqr(p) - iEpsilon) / sqr(q)));
}

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
T divide_finite(T a, T b) {
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

double integrandThreshbnr(double x, double p, double m1, double m2, double q) noexcept
{
  return refnfn(x, p, m1, m2, q);
}

Array<double,1,1> dd(double x, double p, double m1, double m2, double q) noexcept
{
  Array<double,1,1> dydx;
  dydx(0) = -integrandThreshbnr(x, p, m1, m2, q);
  return dydx;
}

// Returns real part of integral
double bIntegral(double p, double m1, double m2, double q)
{
  using namespace flexiblesusy;

  const double from = 0.0, to = 1.0, guess = 0.1, hmin = TOLERANCE * 1.0e-5;
  const double eps = TOLERANCE * 1.0e-3;
  Array<double,1,1> v;
  v(0) = 1.0;

  auto deriv = [p, m1, m2, q] (double x, const Eigen::Array<double,1,1>&) {
     return dd(x, p, m1, m2, q);
  };

  runge_kutta::integrateOdes(v, from, to, eps, guess, hmin,
                             deriv,
                             runge_kutta::odeStepper<Eigen::Array<double,1,1>, decltype(deriv)>);

  return v(0) - 1.0;
}

double fB(const std::complex<double>& a) noexcept
{
  const double x = a.real();

  if (fabs(x) < EPSTOL)
     return -1. - x + sqr(x) * 0.5;

  if (close(x, 1., EPSTOL))
     return -1.;

  return std::real(std::log(1. - a) - 1. - a * std::log(1.0 - 1.0 / a));
}

constexpr double pow3(double a) noexcept { return a*a*a; }
constexpr double pow6(double a) noexcept { return a*a*a*a*a*a; }

} // anonymous namespace

double calcDerivative(double (*func)(double), double x, double h, double
		      *err){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw flexiblesusy::SoftsusyNumericsError("h must be nonzero in numerics.cpp:calcDerivative");


  DoubleMatrix a(NTAB, NTAB);
  hh = h;
  a(1, 1) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
  *err = BIG;
  for (i=2; i<=NTAB; i++) {
    hh /= CON;
    a(1, i) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
    fac = CON2;
    for (j=2; j<=i; j++) {
      a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1.0);
      fac = CON2 * fac;
      errt = maximum(fabs(a(j, i) - a(j-1, i)), fabs(a(j, i) - a(j-1, i-1)));
      if (errt <= *err) {
	*err = errt;
	ans = a(j, i);
      }
    }
    if (fabs(a(i, i) - a(i-1, i-1)) >= SAFE * (*err)) break;
  }

  return ans;
}

double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin)
{
  const double R = 0.61803399, C = 1.0 - R;
  double f1, f2, x0, x1, x2, x3;
  
  x0 = ax; 
  x3 = cx; 
  if (fabs(cx - bx) > fabs(bx - ax)) {
    x1 = bx; 
    x2 = bx + C * (cx - bx); 
  } else {
    x2 = bx; 
    x1 = bx - C * (bx - ax); 
  }
  f1 = (*f)(x1); 
  f2 = (*f)(x2); 
  while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
    if (f2 < f1) {
      shft3(x0, x1, x2, R * x1 + C * x3);
      shft2(f1, f2, (*f)(x2));
    } else {
      shft3(x3, x2, x1, R * x2 + C * x0);
      shft2(f2, f1, (*f)(x1));
	}
  }
  if (f1 < f2) {
    *xmin = x1; 
    return f1; 
  } else {
    *xmin = x2; 
    return f2; 
  }
}

/*
  Analytic expressions follow for above integrals: sometimes useful!
  From hep-ph/9606211
  Note it returns the REAL PART ONLY. 
*/
double b0(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  double b0l = B0(p*p, m1*m1, m2*m2).real();
  //  return B0(p*p, m1*m1, m2*m2).real();
#endif

  // protect against infrared divergence
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
     return 0.0;

  double ans  = 0.;
  const double mMin = minimum(fabs(m1), fabs(m2));
  const double mMax = maximum(fabs(m1), fabs(m2));

  const double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
  /// Try to increase the accuracy of s
  const double dmSq = mMaxSq - mMinSq;
  const double s = pSq + dmSq;

  const double pTest = divide_finite(pSq, mMaxSq);
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-10;

  /// p is not 0  
  if (pTest > pTolerance) {  
     const std::complex<double> iEpsilon(0.0, EPSTOL * mMaxSq);
     const std::complex<double> xPlus =
        (s + sqrt(sqr(s) - 4. * pSq * (mMaxSq - iEpsilon))) / (2. * pSq);
     const std::complex<double> xMinus = 2. * (mMaxSq - iEpsilon) /
        (s + sqrt(sqr(s) - 4. * pSq * (mMaxSq - iEpsilon)));

    ans = -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
  } else {
    if (close(m1, m2, EPSTOL)) {
      ans = - log(sqr(m1 / q));
    } else {
      const double Mmax2 = mMaxSq, Mmin2 = mMinSq;
      if (Mmin2 < 1.e-30) {
	ans = 1.0 - log(Mmax2 / sqr(q));
      } else {
	ans = 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2) 
	  / (Mmin2 - Mmax2);
      }
    }
  }   

#ifdef USE_LOOPTOOLS
  if (!close(b0l, ans, 1.0e-3)) {
    cout << "DEBUG Err: DB0(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b0l/ans << endl;
    cout << "SOFTSUSY  B0=" << ans << endl;
    cout << "LOOPTOOLS B0=" << b0l << endl;
  }
#endif

  return ans;
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
double b1(double p, double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  double b1l = -B1(p*p, m1*m1, m2*m2).real();
  //    return b1l;
#endif

  // protect against infrared divergence
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
     return 0.0;

  double ans = 0.;

  const double p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2), q2 = sqr(q);
  const double pTest = divide_finite(p2, maximum(m12, m22));

  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-4; 

  if (pTest > pTolerance) {
    ans = (a0(m2, q) - a0(m1, q) + (p2 + m12 - m22)
	   * b0(p, m1, m2, q)) / (2.0 * p2);
  } else if (fabs(m1) > 1.0e-15 && fabs(m2) > 1.0e-15) { ///< checked
    const double m14 = sqr(m12), m24 = sqr(m22);
    const double m16 = m12*m14 , m26 = m22*m24;
    const double m18 = sqr(m14), m28 = sqr(m24);
    const double p4 = sqr(p2);
    if (fabs(m12 - m22) < pTolerance * maximum(m12, m22)) {
       ans = 0.08333333333333333*p2/m22
          + 0.008333333333333333*p4/m24
          + sqr(m12 - m22)*(0.041666666666666664/m24 +
                            0.016666666666666666*p2/m26 +
                            0.005357142857142856*p4/m28)
          + (m12 - m22)*(-0.16666666666666666/m22 -
                         0.03333333333333333*p2/m24 -
                         0.007142857142857142*p4/m26)
          - 0.5*log(m22/q2);
    } else {
       ans = (3*m14 - 4*m12*m22 + m24 - 2*m14*log(m12/m22))/(4.*sqr(m12 - m22))
          + (p2*(4*pow3(m12 - m22)*
                 (2*m14 + 5*m12*m22 - m24) +
                 (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
                 - 12*m14*m22*(2*sqr(m12 - m22) + (2*m12 + 3*m22)*p2)*log(m12/m22)))/
          (24.*pow6(m12 - m22)) - 0.5*log(m22/q2);
    }
  } else {
    ans = bIntegral(p, m1, m2, q);
  }

#ifdef USE_LOOPTOOLS
  if (!close(b1l, ans, 1.0e-3)) {
    cout << " Test=" << pTest << " ";
    cout << "DEBUG Err: Db1(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b1l/ans << endl;
    cout << "SOFTSUSY  B1=" << ans << " B0=" << b0(p, m1, m2, q) << endl;
    cout << "LOOPTOOLS B1=" << b1l << " B0=" << B0(p*p, m1*m1, m2*m2).real() 
	 << endl;
  }
#endif
  
  return ans;
}

double b22(double p,  double m1, double m2, double q) {
#ifdef USE_LOOPTOOLS
  setmudim(q*q);
  double b22l = B00(p*p, m1*m1, m2*m2).real();
#endif

  // protect against infrared divergence
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
     return 0.0;

  double answer = 0.;
  
  /// Decides level at which one switches to p=0 limit of calculations
  const double p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2);
  const double pTolerance = 1.0e-10;

  if (p2 < pTolerance * maximum(m12, m22) ) {
    // m1 == m2 with good accuracy
    if (close(m1, m2, EPSTOL)) {
      answer = -m12 * log(sqr(m1 / q)) * 0.5 + m12 * 0.5;
    }
    else
      /// This zero p limit is good
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL) {
	answer = 0.375 * (m12 + m22) - 0.25 *
	  (sqr(m22) * log(sqr(m2 / q)) - sqr(m12) *
	   log(sqr(m1 / q))) / (m22 - m12);
      }
      else
	if (fabs(m1) < EPSTOL) {
	  answer = 0.375 * m22 - 0.25 * m22 * log(sqr(m2 / q));
	}
	else {
	  answer = 0.375 * m12 - 0.25 * m12 * log(sqr(m1 / q));
	}
  }
  else {// checked
    const double b0Save = b0(p, m1, m2, q);

    answer = 1.0 / 6.0 * 
      (0.5 * (a0(m1, q) + a0(m2, q)) + (m12 + m22 - 0.5 * p2)
       * b0Save + (m22 - m12) / (2.0 * p2) *
       (a0(m2, q) - a0(m1, q) - (m22 - m12) * b0Save) +
       m12 + m22 - p2 / 3.0);
  }

#ifdef USE_LOOPTOOLS
  if (!close(b22l, answer, 1.0e-3)) {
    cout << " DEBUG Err: Db22(" << p << ", " << m1 << ", " << m2 
	 << ", "  << q << ")=" << 1.-b22l/answer << endl;
    cout << "SOFTSUSY  B22=" << answer << " B0=" << b0(p, m1, m2, q) << endl;
    cout << "LOOPTOOLS B22=" << b22l << " B0=" << B0(p*p, m1*m1, m2*m2).real() 
	 << endl;
  }
#endif

  return answer;
}

// debugged 23.01.07 - thanks to Shindou Tetsuo 
double d0(double m1, double m2, double m3, double m4) {
  if (close(m1, m2, EPSTOL)) {
    double m2sq = sqr(m2), m3sq = sqr(m3), m4sq = sqr(m4);

    if (close(m2,0.,EPSTOL)) {
       // d0 is undefined for m1 == m2 == 0
       return 0.;
    } else if (close(m3,0.,EPSTOL)) {
       return (-sqr(m2) + sqr(m4) - sqr(m2) * log(sqr(m4/m2)))/
          sqr(m2 * sqr(m2) - m2 * sqr(m4));
    } else if (close(m4,0.,EPSTOL)) {
       return (-sqr(m2) + sqr(m3) - sqr(m2) * log(sqr(m3/m2)))/
          sqr(m2 * sqr(m2) - m2 * sqr(m3));
    } else if (close(m2, m3, EPSTOL) && close(m2, m4, EPSTOL)) {
      return 1.0 / (6.0 * sqr(m2sq));
    } else if (close(m2, m3, EPSTOL)) {
      return (sqr(m2sq) - sqr(m4sq) + 2.0 * m4sq * m2sq * log(m4sq / m2sq)) / 
	(2.0 * m2sq * sqr(m2sq - m4sq) * (m2sq - m4sq));
    } else if (close(m2, m4, EPSTOL)) {
      return (sqr(m2sq) - sqr(m3sq) + 2.0 * m3sq * m2sq * log(m3sq / m2sq)) / 
	(2.0 * m2sq * sqr(m2sq - m3sq) * (m2sq - m3sq));
    } else if (close(m3, m4, EPSTOL)) {
      return -1.0 / sqr(m2sq - m3sq) * 
	((m2sq + m3sq) / (m2sq - m3sq) * log(m3sq / m2sq) + 2.0);
    }

    return 
      (m4sq / sqr(m2sq - m4sq) * log(m4sq / m2sq) + 
       m4sq / (m2sq * (m2sq - m4sq)) -
       m3sq / sqr(m2sq - m3sq) * log(m3sq / m2sq) -
       m3sq / (m2sq * (m2sq - m3sq))) / (m3sq - m4sq);
  }
  return (c0(m1, m3, m4) - c0(m2, m3, m4)) / (sqr(m1) - sqr(m2));
}

double d27(double m1, double m2, double m3, double m4) {// checked

  if (close(m1, m2, EPSTOL)) {
    const double m1n = m1 + TOLERANCE * 0.01;
    return (sqr(m1n) * c0(m1n, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
      / (4.0 * (sqr(m1n) - sqr(m2)));
  }
  return (sqr(m1) * c0(m1, m3, m4) - sqr(m2) * c0(m2, m3, m4)) 
    / (4.0 * (sqr(m1) - sqr(m2)));
}

// Bug-fixed 14.10.02 by T. Watari and collaborators - many thanks!
double c0(double m1, double m2, double m3) {

#ifdef USE_LOOPTOOLS
  double q = 100.;
  setmudim(q*q); 
  double psq = 0.;
  double c0l = C0(psq, psq, psq, m1*m1, m2*m2, m3*m3).real();
#endif

  double ans = 0.;

  if (close(m1,0.,EPSTOL) && close(m2,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m1 == m2 == m3 == 0
     ans = 0.;
  } else if (close(m2,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m2 == m3 == 0
     ans = 0.;
  } else if (close(m1,0.,EPSTOL) && close(m3,0.,EPSTOL)) {
     // c0 is undefined for m1 == m3 == 0
     ans = 0.;
  } else if (close(m1,0.,EPSTOL) && close(m2,0.,EPSTOL)) {
     // c0 is undefined for m1 == m2 == 0
     ans= 0.;
  } else if (close(m1,0.,EPSTOL)) {
     if (close(m2,m3,EPSTOL)) {
        ans = -1./sqr(m2);
     } else {
        ans = (-log(sqr(m2)) + log(sqr(m3)))/(sqr(m2) - sqr(m3));
     }
  } else if (close(m2,0.,EPSTOL)) {
     if (close(m1,m3,EPSTOL)) {
        ans = -1./sqr(m1);
     } else {
        ans = log(sqr(m3/m1))/(sqr(m1) - sqr(m3));
     }
  } else if (close(m3,0.,EPSTOL)) {
     if (close(m1,m2,EPSTOL)) {
        ans = -1./sqr(m1);
     } else {
        ans = log(sqr(m2/m1))/(sqr(m1) - sqr(m2));
     }
  } else if (close(m2, m3, EPSTOL)) {
    if (close(m1, m2, EPSTOL)) {
      ans = ( - 0.5 / sqr(m2) ); // checked 14.10.02
    } else {
      ans = ( sqr(m1) / sqr(sqr(m1)-sqr(m2) ) * log(sqr(m2)/sqr(m1))
               + 1.0 / (sqr(m1) - sqr(m2)) ) ; // checked 14.10.02
    }
  } else if (close(m1, m2, EPSTOL)) {
     ans = ( - ( 1.0 + sqr(m3) / (sqr(m2)-sqr(m3)) * log(sqr(m3)/sqr(m2)) )
             / (sqr(m2)-sqr(m3)) ) ; // checked 14.10.02
  } else if (close(m1, m3, EPSTOL)) {
     ans = ( - (1.0 + sqr(m2) / (sqr(m3)-sqr(m2)) * log(sqr(m2)/sqr(m3)))
             / (sqr(m3)-sqr(m2)) ); // checked 14.10.02
  } else {
     ans = (1.0 / (sqr(m2) - sqr(m3)) *
            (sqr(m2) / (sqr(m1) - sqr(m2)) *
             log(sqr(m2) / sqr(m1)) -
             sqr(m3) / (sqr(m1) - sqr(m3)) *
             log(sqr(m3) / sqr(m1))) );
  }

#ifdef USE_LOOPTOOLS
  if (!close(c0l, ans, 1.0e-3)) {
    cout << " DEBUG Err: C0" << m1 << ", " << m2 
	 << ", "  << m3 << ")=" << 1.-c0l/ans << endl;
    cout << "SOFTSUSY  C0=" << ans << endl;
    cout << "LOOPTOOLS C0=" << c0l << endl;
  }
#endif

  return ans;
}

double dilog(double x) {
  return gm2calc::dilog(x);
}

Complex dilog(const Complex & x) {
  return gm2calc::dilog(x);
}

double fps(double z) {
  if (z < 0.25) {
    double y = sqrt(1.0 - 4.0 * z);
    return 2.0 * z / y * 
      (dilog(1.0 - (1.0 - y) / (2.0 * z)) - 
       dilog(1.0 - (1.0 + y) / (2.0 * z)));
  }
  Complex zz(z);
  Complex y = sqrt(1.0 - 4.0 * zz);

  Complex ans = 2.0 * zz / y * 
    (dilog(1.0 - (1.0 - y) / (2.0 * zz)) - 
     dilog(1.0 - (1.0 + y) / (2.0 * zz)));

  /// answer should always be real
  if (ans.imag() > EPSTOL) throw flexiblesusy::SoftsusyNumericsError("Error in fps");
  return ans.real();
}

double fs(double z) {
  return (2.0 * z - 1) * fps(z) - 2.0 * z * (2.0 + log(z));
}

double ffbar(double z) {
  return z * 0.5 * (2.0 + log(z) - fps(z));
}

DoubleMatrix display3x3RealMixing(double theta12, double theta13, 
				  double theta23, double d) {
  Complex eID(cos(d), sin(d));
  double s12 = sin(theta12);
  double s13 = sin(theta13);
  double s23 = sin(theta23);
  
  double c12 = cos(theta12);
  double c13 = cos(theta13);
  double c23 = cos(theta23);

  DoubleMatrix ckmMatrix(3, 3);
  ckmMatrix(1, 1) = c12 * c13;      
  ckmMatrix(1, 2) = s12 * c13; 

  /// phase factor e^i delta: we'll set it to + or - 1 depending on the sign
  /// of s13
  int pf = 1;
  if (s13 < 0.) pf = -1;
  
  ckmMatrix(1, 3) = pf * s13;
  ckmMatrix(2, 1) = (-s12 * c23 - pf * c12 * s23 * s13);
  ckmMatrix(2, 2) = (c12 * c23 - pf * s12 * s23 * s13);
  ckmMatrix(2, 3) = s23 * c13; 
  ckmMatrix(3, 1) = (s12 * s23 - pf * c12 * c23 * s13); 
  ckmMatrix(3, 2) = (-c12 * s23 - pf * s12 * c23 * s13); 
  ckmMatrix(3, 3) = c23 * c13;

  return ckmMatrix;
  
}

double den(double a, int b) {
  double aa = a;
  int i; for (i=1; i<b; i++) aa = aa * a;
  return 1. / aa;
}

double fin(double mm1, double mm2) {
  if (mm1>mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) +
	    dilog(mm2/mm1) - (mm2*dilog(mm2/mm1))/mm1 -
	    (3.*mm2*log(mm1))/mm1 - log(mm1)*log(mm1 - mm2) +
	    (mm2*log(mm1)*log(mm1 - mm2))/mm1 + (3.*mm2*log(mm2))/mm1 -
	    log(mm1)*log(mm2) + (2.*mm2*log(mm1)*log(mm2))/mm1 +
	    log(mm1 - mm2)*log(mm2) - (mm2*log(mm1 - mm2)*log(mm2))/mm1 - 
	    sqr(PI) * 0.25 +
	    (mm2*sqr(PI))/(12.*mm1) + sqr(log(mm1)) - 
	    (3.*mm2*sqr(log(mm1)))/(2.*mm1) -
	    (mm2*sqr(log(mm2)))/(2.*mm1));
  else if (mm1<mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) -
	    dilog(mm1/mm2) + (mm2*dilog(mm1/mm2))/mm1 -
	    (3*mm2*log(mm1))/mm1 + (3*mm2*log(mm2))/mm1 +
	    (mm2*log(mm1)*log(mm2))/mm1 - log(mm1)*
	    log(-mm1 + mm2)+ (mm2*log(mm1)*log(-mm1 + mm2))/mm1 
	    + log(mm2)*log(-mm1 + mm2) -
	    (mm2*log(mm2)*log(-mm1 + mm2))/mm1 + sqr(PI) / 12. -
	    (mm2*sqr(PI))/(4.*mm1) + sqr(log(mm1)) * 0.5 - 
	    (mm2*sqr(log(mm1)))/mm1 -
	    sqr(log(mm2))/2);
  else return 7.-sqr(PI)/6.;  
}


void rotate(DoubleMatrix & r, DoubleMatrix & qt, int n, int i, float a, 
	    float b) {
  int j;
  double c,fact,s,w,y;
  
  if (a == 0.0) {
    c = 0.0;
    s = (b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    fact = b / a;
    c = sign(1.0/sqrt(1.0+(fact*fact)),a);
    s = fact*c;
  } else {
    fact = a / b;
    s = sign(1.0 / sqrt(1.0 + (fact * fact)), b);
    c = fact * s;
  }
  for (j=i; j<=n; j++) {
    y = r(i, j);
    w = r(i+1, j);
    r(i, j) = c * y - s * w;
    r(i+1, j) = s * y + c * w;
  }
  for (j=1; j<=n; j++) {
    y = qt(i, j);
    w = qt(i+1, j);
    qt(i, j) = c * y - s * w;
    qt(i+1, j) = s * y + c * w;
  }
}
