
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include <Eigen/Dense>
#include "numerics.h"
#include "dilog.hpp"
#include "rk.hpp"
#ifdef USE_LOOPTOOLS
#include "clooptools.h"
#endif
#include <boost/bind.hpp>

using namespace softsusy;
using namespace Eigen;

double lnLPoisson(unsigned k, double lambda) {
  if (lambda < 6.66e-66) throw("In lnPoisson: lambda not positive\n");
  double lnL = -lambda + double(k) * log(lambda);
  for (unsigned i=2; i<=k; i++) lnL -= log(double(i));

  return lnL;
}

double LPoisson(unsigned k, double lambda) {
  if (lambda < 6.66e-66) throw("In lnPoisson: lambda not positive\n");
  double L = exp(-lambda) * pow(lambda,int(k));
  for (unsigned i=2; i<=k; i++) L /= double(i);

  return L;
}

double calcDerivative(double (*func)(double), double x, double h, double
		      *err){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw "h must be nonzero in numerics.cpp:calcDerivative";


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

inline void shft2(double & a, double & b, double c) { a = b; b = c; }

inline void shft3(double & a, double & b, double & c, double d) { 
  a = b; b = c; c = d;
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



ArrayXd dd(double x, const ArrayXd& /* y */) {
  ArrayXd dydx(1);
  dydx(0) = -integrandThreshbnr(x);
  return dydx;
}

ArrayXd dd_threadsave(double x, const ArrayXd&, int n1, double p, double m1, double m2, double mt)
{
  ArrayXd dydx(1);
  dydx(0) = -integrandThreshbnr(x, n1, p, m1, m2, mt);
  return dydx;
}

double integrandThreshbnr(double x) {
  return fnfn(x).real();
}

double integrandThreshbnr(double x, int n1, double p, double m1, double m2, double mt)
{
  return fnfn(x, n1, p, m1, m2, mt).real();
}

// Integration routine needs these variables
static double mtInt, pInt, m1Int, m2Int;
static int nInt;

Complex fnfn(double x) {
  const static Complex iEpsilon(0.0, TOLERANCE * 1.0e-20);
  
  double xn = 1.0;
  int i; for (i=1; i<=nInt; i++) xn = xn * x;
  return xn * 
    log( ((1 - x) * sqr(m1Int) + x * sqr(m2Int) - x * (1 - x) *
	  sqr(pInt) - iEpsilon)
	 / sqr(mtInt));
}

Complex fnfn(double x, int n1, double p, double m1, double m2, double mt)
{
  const static Complex iEpsilon(0.0, TOLERANCE * 1.0e-20);

  double xn = 1.0;

  for (int i = 1; i <= n1; i++)
     xn *= x;

  return xn *
    log(((1 - x) * sqr(m1) + x * sqr(m2)
         - x * (1 - x) * sqr(p) - iEpsilon) / sqr(mt));
}

DoubleVector dilogarg(double t, const DoubleVector & /* y */) {

  const double eps = TOLERANCE * 1.0e-20;

  DoubleVector dydx(1);
  dydx(1) = -log(fabs(1 - t + eps)) / (t + eps);

  return dydx;
}

namespace {

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
T divide_finite(T a, T b) {
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

} // anonymous namespace

/*
double dilog(double x) {
  // Set global variables so that integration function can access them
  double from = 0.0, to = x, guess = 0.1, hmin = TOLERANCE * 1.0e-5;

  DoubleVector v(1); 
  double eps = TOLERANCE * 1.0e-5;
  v(1) = 1.0; 

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dilogarg, odeStepper); 
  
  return v(1) - 1.0;
}
*/

// Returns real part of integral
double bIntegral(int n1, double p, double m1, double m2, double mt) {
  using namespace flexiblesusy::runge_kutta;

  // Set global variables so that integration function can access them
  nInt = n1; pInt = p; m1Int = m1; m2Int = m2; mtInt = mt;
  double from = 0.0, to = 1.0, guess = 0.1, hmin = TOLERANCE * 1.0e-5;
  
  ArrayXd v(1); double eps = TOLERANCE * 1.0e-3;
  v(0) = 1.0;

  // Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  // odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, dd, odeStepper);
  
  return v(0) - 1.0;
}

// Returns real part of integral
double bIntegral_threadsave(int n1, double p, double m1, double m2, double mt) {
  using namespace flexiblesusy;

  const double from = 0.0, to = 1.0, guess = 0.1, hmin = TOLERANCE * 1.0e-5;
  const double eps = TOLERANCE * 1.0e-3;
  ArrayXd v(1);
  v(0) = 1.0;

  runge_kutta::Derivs derivs = boost::bind(&dd_threadsave, _1, _2, n1, p, m1, m2, mt);

  runge_kutta::integrateOdes(v, from, to, eps, guess, hmin, derivs,
                             runge_kutta::odeStepper);

  return v(0) - 1.0;
}

double fB(const Complex & a) {
  /// First, special cases at problematic points
  const double x = a.real();
  if (fabs(x) < EPSTOL) {
    const double ans = -1. - x + sqr(x) * 0.5;
    return ans;
  }
  if (close(x, 1., EPSTOL)) return -1.;

  const Complex ans(log(1. - a) - 1. - a * log(1.0 - 1.0 / a));

  return ans.real();
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
    const Complex iEpsilon(0.0, EPSTOL * sqr(mMax));
    
    Complex xPlus, xMinus;

    xPlus = (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon))) /
      (2. * sqr(p));
    xMinus = 2. * (sqr(mMax) - iEpsilon) / 
      (s + sqrt(sqr(s) - 4. * sqr(p) * (sqr(mMax) - iEpsilon)));

    ans = -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
  } else {
    if (close(m1, m2, EPSTOL)) {
      ans = - log(sqr(m1 / q));
    } else {
      const double Mmax2 = sqr(mMax), Mmin2 = sqr(mMin);
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

double fB_fast(const Complex& a) {

  const double x = a.real();

  if (fabs(x) < EPSTOL) {
    return -1. - x + sqr(x) * 0.5;
  }

  if (close(x, 1., EPSTOL))
     return -1.;

  return Complex(log(1. - a) - 1. - a * log(1.0 - 1.0 / a)).real();
}

double b0_fast(double p, double m1, double m2, double q) {
  // protect against infrared divergence
  if (close(p, 0.0, EPSTOL) && close(m1, 0.0, EPSTOL)
      && close(m2, 0.0, EPSTOL))
     return 0.0;

  double ans  = 0.;
  const double mMin = minimum(fabs(m1), fabs(m2));
  const double mMax = maximum(fabs(m1), fabs(m2));

  const double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax),
     q2 = sqr(q);
  /// Try to increase the accuracy of s
  const double dmSq = mMaxSq - mMinSq;
  const double s = pSq + dmSq, s2 = sqr(s);

  const double pTest = divide_finite(pSq, mMaxSq);
  /// Decides level at which one switches to p=0 limit of calculations
  const double pTolerance = 1.0e-6;

  /// p is not 0
  if (pTest > pTolerance) {
     const Complex iEpsilon(0.0, EPSTOL * mMaxSq);
     if (mMinSq < 1.e-30) {
        if (mMaxSq < 1.e-30) {
           ans = - log((-pSq - iEpsilon)/q2).real() + 2.;
        } else if (close(pSq, mMaxSq, EPSTOL)) {
           ans = - log(mMaxSq / q2) + 2.;
        } else {
           ans = - log(mMaxSq / q2) + 2. + (mMaxSq - pSq) / pSq * log((mMaxSq - pSq - iEpsilon)/mMaxSq).real();
        }
     } else if (close(pSq, mMaxSq, EPSTOL) && close(pSq, mMinSq, EPSTOL)) {
        ans = - log(mMaxSq / q2) + 2. - PI/sqrt(3.);
     } else {
        const Complex xPlus(0.5 * (s + sqrt(s2 - 4. * pSq * (mMaxSq - iEpsilon))) / pSq);
        const Complex xMinus(2. * (mMaxSq - iEpsilon) / (s + sqrt(s2 - 4. * pSq * (mMaxSq - iEpsilon))));
        ans = -2.0 * log(p / q) - fB_fast(xPlus) - fB_fast(xMinus);
     }
  } else {
     if (close(m1, m2, EPSTOL)) {
        ans = - log(sqr(m1 / q));
     } else {
        const double Mmax2 = mMaxSq, Mmin2 = mMinSq;
        if (Mmin2 < 1.e-30) {
           ans = 1.0 - log(Mmax2 / q2);
        } else {
           ans = 1.0 - log(Mmax2 / q2) + Mmin2 * log(Mmax2 / Mmin2)
              / (Mmin2 - Mmax2);
        }
     }
  }

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
  } else if (fabs(m1) > 1.0e-15 && !close(m1, m2, EPSTOL)
	     && fabs(m2) > 1.0e-15) { ///< checked
    ans = 0.5 * (1. + log(q2 / m22) +
		 sqr(m12 / (m12 - m22)) * log(m22 / m12) +
		 0.5 * (m12 + m22) / (m12 - m22)
		 );
  } else {
    ans = bIntegral_threadsave(1, p, m1, m2, q);
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

    if (close(m2, m3, EPSTOL) && close(m2, m4, EPSTOL)) 
      return 1.0 / (6.0 * sqr(m2sq));
    else if (close(m2, m3, EPSTOL)) {
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

  double ans;

  if (close(m2, m3, EPSTOL)) {
    if (close(m1, m2, EPSTOL)) {
      ans = ( - 0.5 / sqr(m2) ); // checked 14.10.02
    }
    else {
      ans = ( sqr(m1) / sqr(sqr(m1)-sqr(m2) ) * log(sqr(m2)/sqr(m1))
               + 1.0 / (sqr(m1) - sqr(m2)) ) ; // checked 14.10.02
    }
  }
  else
    if (close(m1, m2, EPSTOL)) {
      ans = ( - ( 1.0 + sqr(m3) / (sqr(m2)-sqr(m3)) * log(sqr(m3)/sqr(m2)) )
               / (sqr(m2)-sqr(m3)) ) ; // checked 14.10.02
    }
    else
      if (close(m1, m3, EPSTOL)) {
        ans = ( - (1.0 + sqr(m2) / (sqr(m3)-sqr(m2)) * log(sqr(m2)/sqr(m3))) 
                 / (sqr(m3)-sqr(m2)) ); // checked 14.10.02
      }
      else {
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

double truncGaussWidthHalf(long & idum) {
  for (;;) {
    double a = ran1(idum);
    if (a < 0.5) return a * 2.0;
	else return fabs(gasdev(idum)) + 1.;
  }
}

double gasdev(long & idum) {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

double ran1(long & idum) {
  const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32, 
    NDIV = 1+(IM-1)/NTAB;
  const double AM = 1.0 / double(IM), EPS = 1.2e-15, RNMX = 1.0 - EPS;
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (idum <= 0 || !iy) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(idum)/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

double cauchyRan(long & idum) {
  double x = ran1(idum) - 0.5;
  double unNormalised = tan(x * PI);
  return unNormalised;
}

int bin(double data, double start, double end, int numBins) {
  double range = end - start;
  double binSize = range / double(numBins);
  return int((data - start) / binSize + 1.);
}

double logOfSum(double a, double b) {
  double max = maximum(a, b);
  double min = minimum(a, b);

  if (max + min < 0. || max == 0.) return asin(1.0);
  double ans = log(max);
  ans = ans + log (1.0 + min / max);

  return ans;
}

double sumOfExp(double a, double b) {
  double max = maximum(a, b);
  double min = minimum(a, b);
  
  /// this underflow determines whether it's actually worth adding the two
  /// figures... any more than 15 digits and it's just going to be the maximum
  /// one. 
  const double underflowTol = 15.0 * log(10.);
  if ((max - min) > underflowTol) return exp(max);
  else return exp(max) * (1.0 + exp(min - max));
}

double dilog(double x) {
  return gm2calc::dilog(x);
}

Complex dilog(const Complex & x) {
  return gm2calc::dilog(x);
}


DoubleVector getRandomDirection(int n, int & numChanged,
				long & idum) {
  DoubleVector del(n);

  int i; 

  double length = 0.;
  do {
    numChanged = 0;
    for (i=1; i<=n; i++) {
	del(i) = gasdev(idum);
	numChanged++;
    }
    length = sqrt(del.dot(del));
  }
  while (length < EPSTOL);

  del = 1.0 / sqrt(del.dot(del)) * del;
  return del;
}


double calcCL(double cl, const DoubleVector & l) {
  DoubleVector s(l.sort());

  double fr = 0., likelihood = 0., tot = 0.;
  int i; 
  for (i=s.displayStart(); i<=s.displayEnd(); i++) tot = tot + s.display(i);
  for (i=s.displayStart(); i<=s.displayEnd(); i++) {
    fr = fr + s.display(i) / tot;
    if (fr > (1.0 - cl)) {
      likelihood = s.display(i);
      // bail out of loop
      i = s.displayEnd() + 1;
    }
  }
  return likelihood;
}

int *ivector(long nl, long nh) {
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
	if (!v) throw("allocation failure in ivector()\n");
	return v-nl+1;
}

/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long /* nh */) {
	free((char *) (v+nl-1));
}

double calc1dFraction(double y, const DoubleVector & l) {
  double n = 0., lMax = l.max();
  double tot = 0.;
  int i; for (i=1; i<=l.displayEnd(); i++) {
    if (l.display(i) <= y * lMax) n = n + l.display(i);
    tot = tot + l.display(i);
  }

  return n / tot;
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
  if (ans.imag() > EPSTOL) throw("Error in fps\n");
  return ans.real();
}

double fs(double z) {
  return (2.0 * z - 1) * fps(z) - 2.0 * z * (2.0 + log(z));
}

double ffbar(double z) {
  return z * 0.5 * (2.0 + log(z) - fps(z));
}

#define JMAX 20

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n, double /* EPS */) {
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

double qtrap(double (*func)(double), double a, double b, double EPS) {
  //	double midpnt(double (*func)(double), double a, double b, int n);
  //	void nrerror(char error_text[]);

	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=midpnt(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		if (s == 0.0 && olds == 0.0 && j > 6) return s;
		olds=s;
	}
	throw("Too many steps in routine qtrap\n");
	return 0.0;
}
#undef JMAX

double midpnt(double (*func)(double), double a, double b, int n) {
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double edgefn(double x, double y, double z) {
  return (x - y) * (y - z) / y;
}
double mllMax(double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double expression = edgefn(xi, l, chi);
  if (mChi2 < mSlep) expression = sqr(mChi2-mChi1);
  if (expression < 0.) return -numberOfTheBeast; 
  return sqrt(expression);
}

double mllq(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  if (sqr(l) < q * chi && q * chi < sqr(xi) && sqr(xi) * chi < q * sqr(l))
    return (mSq - mChi1);
  double exp1  = edgefn(q, xi, chi);
  if (exp1 < 0.) exp1 = -numberOfTheBeast; 
  double exp2  = edgefn(q, l, chi);
  if (exp2 < 0.) exp2 = -numberOfTheBeast; 
  double exp3  = (q * l - xi * chi) * (xi - l) / (xi * l);
  if (exp3 < 0.) exp3 = -numberOfTheBeast; 
  double expression = maximum(maximum(exp1, exp2), exp3);
  if (mChi2 < mSlep  && xi > sqr(q) * sqr(chi) ) return (mSq - mChi1);
  if (mChi2 < mSlep  && xi < sqr(q) * sqr(chi) )  expression = edgefn(q,xi,chi);  
  return sqrt(expression);
}

double llqThresh(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double ans   = 2.0 * l * (q - xi) * (xi - chi) + 
    (q + xi) * (xi - l) * (l - chi) - (q - xi) * 
    sqrt(sqr(xi + l) * sqr(l + chi) - 16. * xi * sqr(l) * chi); 
  if (ans < 0.) return -numberOfTheBeast;
  return sqrt(ans / (4.0 * l * xi));
}

double lqnear(double mSq, double /* mChi1 */, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double q     = sqr(mSq);
  double expression = (q - xi) * (xi - l) / xi;
  if (expression < 0.) return -numberOfTheBeast;
  return sqrt(expression);
}

double lqfar(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double expression = (q - xi) * (l - chi) / l;
  if (expression < 0.) return -numberOfTheBeast;
  return sqrt(expression);
}

double lqhigh(double mSq, double mChi1, double mSlep, double mChi2) {
  double mlqnear = lqnear(mSq, mChi1, mSlep, mChi2);
  double mlqfar  = lqfar(mSq, mChi1, mSlep, mChi2);
if (mChi2 < mSlep) {return sqrt(edgefn(sqr(mSq),sqr(mChi2),sqr(mChi1)));}
  return maximum(mlqnear, mlqfar);
}

double lqlow(double mSq, double mChi1, double mSlep, double mChi2) {
  double xi    = sqr(mChi2);
  double l     = sqr(mSlep);
  double chi   = sqr(mChi1);
  double q     = sqr(mSq);
  double mlqnear = lqnear(mSq, mChi1, mSlep, mChi2);
  double mlqother  = (q - xi) * (l - chi) / (2.0 * l - chi);
  if (mlqother < 0.) mlqother = -numberOfTheBeast;
  else mlqother = sqrt(mlqother);

  return minimum(mlqnear, mlqother);
}



#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b, double EPS) {
  double ss = 0.,dss = 0.;
  DoubleVector s(JMAXP+1),h(JMAXP+1);
  int j;
  
  h(1) = 1.0;
  for (j=1; j<=JMAX; j++) {
    s(j)=trapzd(func,a,b,j);
    if (j >= K) {
      DoubleVector h2(K), s2(K);
      for (unsigned int ii=1; ii<=K; ii++) {
	h2(ii) = h(j-K+ii); 
	s2(ii) = s(j-K+ii);
      }	
      polint(h2,s2,0.0, ss, dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    s(j+1)=s(j);
    h(j+1)=0.25 * h(j);
  }
  throw("Too many steps in routine qromb\n");
  return 0.0;
}

double qromb2(double (*func)(double), double a, double b, double EPS) {
  double ss = 0.,dss = 0.;
  DoubleVector s(JMAXP+1),h(JMAXP+1);
  int j;
  
  h(1) = 1.0;
  for (j=1; j<=JMAX; j++) {
    s(j)=trapzd(func,a,b,j);
    if (j >= K) {
      DoubleVector h2(K), s2(K);
      for (unsigned int ii=1; ii<=K; ii++) {
	h2(ii) = h(j-K+ii); 
	s2(ii) = s(j-K+ii);
      }	
      polint(h2,s2,0.0, ss, dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    s(j+1)=s(j);
    h(j+1)=0.25 * h(j);
  }
  throw("Too many steps in routine qromb\n");
  return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K

void polint(const DoubleVector&  xa, const DoubleVector & ya, double x, 
	    double & y, double & dy) {
  if (xa.displayEnd() != ya.displayEnd()) throw("Incompatible vector lengths in numerics.cpp:polint\n");
  
  int n = xa.displayEnd();
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  DoubleVector c(n), d(n);

  dif=fabs(x-xa.display(1));
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa.display(i))) < dif) {
      ns=i;
      dif=dift;
    }
    c(i)=ya.display(i);
    d(i)=ya.display(i);
  }
  y=ya.display(ns--);
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa.display(i)-x;
      hp=xa.display(i+m)-x;
      w=c(i+1)-d(i);
      if ( (den=ho-hp) == 0.0) throw("Error in routine polint\n");
      den=w/den;
      d(i)=hp*den;
      c(i)=ho*den;
    }
    y += (dy=(2*ns < (n-m) ? c(ns+1) : d(ns--)));
  }
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

void getAngles(const DoubleMatrix & v, double & t12, double & t13, 
	       double & t23, double & d) {
  int pf = 1;
  t23 = atan(v.display(2, 3) / v.display(3, 3));
  t13 = asin(v.display(1, 3));
  if (t13 < 0.) {
    pf = -1; 
    t13 = -t13;
  }
  d = acos(double(pf));
  t12 = asin(v.display(1, 2) / cos(t13));
}

bool midPtStep(DoubleVector & xi, 
	      DoubleVector (*derivs)(double t, const DoubleVector & v), 
	      double tInitial, double tStep) {
  /// initial guess
  DoubleVector xiPlus1Old(xi + tStep * derivs(tInitial + tStep, xi));
  DoubleVector xiPlus1New(xi.displayEnd());

  double diff = numberOfTheBeast, delta = 0.;
  int count = 0;
  while (diff > 1.0e-15 && count < 100) {
    count++;
    xiPlus1New = xi + tStep * derivs(tInitial + tStep, 0.5*(xi + xiPlus1Old));
    double max = -1.; /// difference between iterations
      int i; for (i=1; i<=xi.displayEnd(); i++) {
	if (fabs(xiPlus1New(i)) < 1.0e-12) 
	delta = fabs(xiPlus1New(i) - xiPlus1Old(i));
      else delta = fabs(1.- xiPlus1Old(i) / xiPlus1New(i));
      if (delta > max) max = delta;
      }

    diff = max;

    xiPlus1Old = xiPlus1New;
  }
  xi = xiPlus1New;
  if (count < 100) return false;
  return true;
}

/// Do a fixed number of steps of approximately reversible integration: hoping
/// it will help convergence in difficult cases. It is really *slow* though.
bool integrateReversibly(DoubleVector & xi, 
			 DoubleVector (*derivs)(double t, 
						const DoubleVector & v), 
			 double tInitial, double tFinal, int numSteps) {
  bool err = false;
  double tStep = (tFinal - tInitial) / double(numSteps);

  double t = tInitial;
  int i; for (i=0; i<=numSteps; i++) {

    if (midPtStep(xi, derivs, t, tStep)) err = true;
    t = tStep * double(i) + tInitial;
  }
  return err;
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

double zriddr(double (*func)(double), double x1, double x2, double xacc) {
  const int MAXIT = 60;
  const double UNUSED = -1.11e30;
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  
  fl=(*func)(x1);
  fh=(*func)(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=UNUSED;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm=(*func)(xm);
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0) return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) return ans;
      ans=xnew;
      fnew=(*func)(ans);
      if (fnew == 0.0) return ans;
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else throw("In zriddr: never get here.");
      if (fabs(xh-xl) <= xacc) return ans;
    }
    throw("zriddr exceeded maximum iterations\n");
  }
  else {
    if (close(fl, 0.0, EPSTOL)) return x1;
    if (close(fh, 0.0, EPSTOL)) return x2;
    ostringstream s;
    s << "Root must be bracketed in zriddr. f(" << x1 << ")=" << fl 
      << ", fh(" << x2 << ")=" << fh << endl;
    throw(s.str());
  }
  return 0.0;
}

/// You will need to clear this lot up....
DoubleMatrix fdjac(int n, DoubleVector x, const DoubleVector & fvec,
	   void (*vecfunc)(const DoubleVector &, DoubleVector &)) {
  double EPS = maximum(TOLERANCE, 1.0e-4);
  int i,j;
  double h,temp;
  
  DoubleVector f(1, n);
  DoubleMatrix df(n, n);
  for (j=1; j<=n; j++) {
    temp = x(j);
    h = EPS * fabs(temp);
    if (h == 0.0) h = EPS;
    x(j) = temp + h;
    h = x(j) - temp;
    (*vecfunc)(x, f);
    x(j) = temp;
    for (i=1; i<=n; i++) df(i, j) = (f(i) - fvec.display(i)) / h;
  }
  return df;
}

bool lnsrch(const DoubleVector & xold, double fold, const DoubleVector & g, 
	    DoubleVector & p, 
	    DoubleVector & x, double & f, double stpmax, 
	    void (*vecfunc)(const DoubleVector &, DoubleVector &), 
	    DoubleVector & fvec) {
  double ALF = TOLERANCE;
  double TOLX = TOLERANCE * 1.0e-3;
  
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0,fold2=0.0,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
  
  bool err = false;
  sum = sqrt(p.dot(p));
  if (sum > stpmax) p = (stpmax / sum) * p;

  for (slope=0.0, i=1; i<=xold.displayEnd(); i++)
    slope += g(i) * p(i);
  test = 0.0;
  for (i=1; i<=xold.displayEnd(); i++) {
    temp=fabs(p(i)) / maximum(fabs(xold(i)), 1.0);
    if (temp > test) test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;) {
    x = xold + alam * p;
    vecfunc(x, fvec); 
    f = fvec.dot(fvec);
    if (alam < alamin) {
      x = xold;
      err = true;
      return err;
    } else if (f <= fold + ALF * alam * slope) return err;
    else {
      if (alam == 1.0)
	tmplam = -slope / (2.0 * (f - fold - slope));
      else {
	rhs1 = f - fold - alam * slope;
	rhs2 = f2 - fold2 - alam2 * slope;
	a=(rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
	b=(-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / 
	  (alam - alam2);
	if (a == 0.0) tmplam = -slope / (2.0 * b);
	else {
	  disc = b * b - 3.0 * a * slope;
	  if (disc < 0.0) throw("Roundoff problem in lnsrch.\n");
	  else tmplam = (-b + sqrt(disc)) / (3.0 * a);
	}
	if (tmplam > 0.5 * alam)
	  tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = f;
    fold2 = fold;
    alam = maximum(tmplam, 0.1 * alam);
  }
  return err;
}

/// Get rid of int n
void lubksb(const DoubleMatrix & a, int n, int *indx, DoubleVector & b) {
  int i,ii=0,ip,j;
  float sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b(ip);
    b(ip)=b(i);
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a(i, j) * b(j);
    else if (sum) ii=i;
    b(i)=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b(i);
    for (j=i+1;j<=n;j++) sum -= a(i, j)*b(j);
    b(i)=sum/a(i, i);
  }
}


/// Get rid of int n
void ludcmp(DoubleMatrix & a, int n, int *indx, double & d) {
  const double TINY = 1.0e-20;
  int i,imax=0,j,k;
  float big,dum,sum=0.0,temp;
  DoubleVector vv(n);

  d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a(i, j))) > big) big=temp;
    if (big == 0.0) {
      ostringstream ii; ii.setf(ios::scientific, ios::floatfield);
      ii.precision(6);
      ii << "Singular matrix in routine ludcmp: ";
      ii << a;
      throw (ii.str());
    }
    vv(i)=1.0 / big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a(i, j);
      for (k=1;k<i;k++) sum -= a(i, k)*a(k, j);
      a(i, j)=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a(i, j);
      for (k=1;k<j;k++)
	sum -= a(i, k)*a(k, j);
      a(i, j)=sum;
      if ( (dum=vv(i)*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a(imax, k);
	a(imax, k)=a(j, k);
	a(j, k)=dum;
      }
      d = -(d);
      vv(imax)=vv(j);
    }
    indx[j]=imax;
    if (a(j, j) == 0.0) a(j, j) = TINY;
    if (j != n) {
      dum=1.0/(a(j, j));
      for (i=j+1;i<=n;i++) a(i, j) *= dum;
    }
  }
}


/// More work can be done on this: get rid of int n and in subfunctions too
bool newt(DoubleVector & x, 
	  void (*vecfunc)(const DoubleVector &, DoubleVector &)) {
  bool err = false; 
  const int MAXITS = 200;    ///< max iterations
  double TOLF   = TOLERANCE; ///< convergence on function values
  double TOLMIN = TOLF * 1.e-2; ///< spurious convergence to min of fmin
  double TOLX   = TOLF * 1.e-3; ///< maximum dx convergence criterion
  const double STPMX  = 100.0; ///< maximum step length allowed in line searches
   
  int n = x.displayEnd();
  int i,its,j,*indx;
  double d,den,f,fold,stpmax,sum=0.0,temp,test;
  
  indx = ivector(1, n);
  DoubleMatrix fjac(n, n);
  DoubleVector g(n), p(n), xold(n);
  DoubleVector fvec(n);

  vecfunc(x, fvec); 
  f = 0.5 * fvec.dot(fvec);
  test = 0.0;
  test = fvec.apply(fabs).max();
  if (test < 0.01 * TOLF) {
    err = false;
    free_ivector(indx, 1, n); return err;
  }
  sum += x.dot(x);
  stpmax = STPMX * maximum(sqrt(sum), (double) n);
  for (its=1; its<=MAXITS; its++) {
    fjac = fdjac(n, x, fvec, vecfunc);
    for (i=1;i<=n;i++) {
      for (sum=0.0, j=1; j<=n; j++) sum += fjac(j, i) * fvec(j);
      g(i) = sum;
    }
    for (i=1; i<=n; i++) xold(i) = x(i);
    fold = f;
    for (i=1; i<=n; i++) p(i) = -fvec(i);
    ludcmp(fjac, n, indx, d);
    lubksb(fjac, n, indx, p);
    err = lnsrch(xold, fold, g, p, x, f, stpmax, vecfunc, fvec);
    test = 0.0;
    for (i=1; i<=n; i++)
      if (fabs(fvec(i)) > test) test = fabs(fvec(i));
    if (test < TOLF) {
      err = false;
      free_ivector(indx, 1, n);
      return err;
    }
    if (err) {
      test = 0.0;
      den = maximum(f, 0.5 * n);
      for (i=1; i<=n; i++) {
	temp = fabs(g(i)) * maximum(fabs(x(i)), 1.0) / den;
	if (temp > test) test = temp;
      }
      err = (test < TOLMIN ? true : false);
      free_ivector(indx, 1, n);
      return err;
    }
    /// Check points aren't getting too close together
    test = 0.0;
    for (i=1; i<=n; i++) {
      temp = (fabs(x(i) - xold(i))) / maximum(fabs(x(i)), 1.0);
      if (temp > test) test = temp;
    }
    if (test < TOLX) { 
      free_ivector(indx, 1, n); 
      return err; 
    }
  }

  return true;
}

DoubleVector testDerivs(double /* x */, const DoubleVector & y) {
  DoubleVector dydx(3);
  dydx(1) = y(1) * y(2) * y(2);
  dydx(2) = y(2) * y(1) * y(3);
  dydx(3) = y(1);
  return dydx;
}

void broydn(DoubleVector x, int & check, 
	    void (*vecfunc)(const DoubleVector &, DoubleVector &)) {
  const int MAXITS = 200;
  double TOLF =  TOLERANCE;
  const double EPS = TOLF * 1.0e-3;
  const double TOLX = EPS;
  const double STPMX = 100.0;
  double TOLMIN = TOLERANCE * 1.0e-2;

  int n = x.displayEnd();
 
  /*	void fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float fmin(float x[]);
	void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void qrdcmp(float **a, int n, float *c, float *d, int *sing);
	void qrupdt(float **r, float **qt, int n, float u[], float v[]);
	void rsolv(float **a, int n, float d[], float b[]);
	int i,its,j,k,restrt,sing,skip;
	float den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
	float *g,*p,**qt,**r,*s,*t,*w,*xold; */

  DoubleVector c(n), d(n), fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n), 
    fvec(n);
  DoubleMatrix qt(n, n), r(n, n);

  vecfunc(x, fvec); 
  double f = 0.5 * fvec.dot(fvec);
  
  double test = fvec.apply(fabs).max();
  if (test < 0.01 * TOLF) {
    check = 0;
    return;
  }
  double sum = x.dot(x);
  double stpmax = STPMX * maximum(sqrt(sum), (double)n);
  int restrt = 1, sing;
  for (int its=1; its<=MAXITS; its++) {
    if (restrt) {
      r = fdjac(n, x, fvec, vecfunc);
      qrdcmp(r, n, c, d, sing);
      if (sing) throw("singular Jacobian in broydn\n");
      for (int i=1; i<=n; i++) {
	for (int j=1; j<=n; j++) qt(i, j) = 0.0;
	qt(i, i) = 1.0;
      }
      for (int k=1; k<n; k++) {
	if (c(k)) {
	  for (int j=1; j<=n; j++) {
	    sum = 0.0;
	    for (int i=k; i<=n; i++)
	      sum += r(i, k) * qt(i, j);
	    sum /= c(k);
	    for (int i=k; i<=n; i++)
	      qt(i, j) -= sum*r(i, k);
	  }
	}
      }
      for (int i=1; i<=n; i++) {
	r(i, i) = d(i);
	for (int j=1; j<i; j++) r(i, j) = 0.0;
      }
    } else {
      s = x - xold;
      int j;
      for (int i=1; i<=n; i++) {
	for (sum=0.0, j=i; j<=n; j++) sum += r(i, j) * s(j);
	t(i) = sum;
      }
      int skip = 1;
      for (int i=1; i<=n; i++) {
	for (sum=0.0, j=1; j<=n; j++) sum += qt(j, i) * t(j);
	w(i) = fvec(i) - fvcold(i) - sum;
	if (fabs(w(i)) >= EPS * (fabs(fvec(i))+fabs(fvcold(i)))) skip=0;
	else w(i)=0.0;
      }
      if (!skip) {
	for (int i=1; i<=n; i++) {
	  for (sum=0.0, j=1; j<=n; j++) sum += qt(i, j) * w(j);
	  t(i) = sum;
	}
	double den = 0.; den += s.dot(s);
	s = s * (1. / den);
	qrupdt(r,qt,n,t,s);
	for (int i=1;i<=n;i++) {
	  if (r(i, i) == 0.0) throw("r singular in broydn\n");
	  d(i) = r(i, i);
	}
      }
    }
    int j;
    for (int i=1; i<=n; i++) {
      for (sum=0., j=1 ; j<=n; j++) sum += qt(i, j) * fvec(j);
      g(i) = sum;
    }
    for (int i=n; i>=1; i--) {
      for (sum=0.0, j=1; j<=i; j++) sum += r(j, i) * g(j);
      g(i) = sum;
    }
    xold = x; 
    fvcold = fvec;
    double fold = f;
    for (int i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt(i, j) * fvec(j);
      p(i) = -sum;
    }
    rsolv(r, n, d, p); 
    DoubleVector fvec(n);
    check = lnsrch(xold, fold, g, p, x, f, stpmax, vecfunc, fvec);
    test = 0.0;
    for (int i=1; i<=n; i++)
      if (fabs(fvec(i)) > test) test = fabs(fvec(i));
    if (test < TOLF) {
      check=0;
      return;
    }
    if (check) {
      if (restrt) return;
      else {
	test = 0.0;
	double den = maximum(f, 0.5 * n);
	for (int i=1;i<=n;i++) {
	  double temp = fabs(g(i)) * maximum(fabs(x(i)), 1.0) / den;
	  if (temp > test) test = temp;
	}
	if (test < TOLMIN) return;
	else restrt = 1;
      }
    } else {
      restrt = 0;
      test = 0.0;
      for (int i=1; i<=n; i++) {
	double temp = (fabs(x(i) - xold(i))) / maximum(fabs(x(i)), 1.0);
	if (temp > test) test = temp;
      }
      if (test < TOLX) return;
    }
  }
  throw("MAXITS exceeded in broydn\n");
}

void qrdcmp(DoubleMatrix & a, int n, DoubleVector & c, DoubleVector & d, 
	    int & sing) {
  int i,j,k;
  double scale,sigma,sum,tau;
  
  sing = 0;
  for (k=1; k<n; k++) {
    scale = 0.0;
    for (i=k; i<=n; i++) scale = maximum(scale, fabs(a(i, k)));
    if (scale == 0.0) {
      sing = 1;
      c(k) = d(k) = 0.0;
    } else {
      for (i=k;i<=n;i++) a(i, k) /= scale;
      for (sum=0.0,i=k;i<=n;i++) sum += sqr(a(i, k));
      sigma = sign(sqrt(sum), a(k, k));
      a(k, k) += sigma;
      c(k) = sigma * a(k, k);
      d(k) = -scale*sigma;
      for (j=k+1; j<=n; j++) {
	for (sum=0.0,i=k;i<=n;i++) sum += a(i, k) * a(i, j);
	tau = sum / c(k);
	for (i=k; i<=n; i++) a(i, j) -= tau * a(i, k);
      }
    }
  }
  d(n) = a(n, n);
  if (d(n) == 0.0) sing = 1;
}


void qrupdt(DoubleMatrix & r, DoubleMatrix & qt, int n, 
	    DoubleVector & u, DoubleVector & v) {
  //	void rotate(float **r, float **qt, int n, int i, float a, float b);
	int i,j,k;

	for (k=n; k>=1; k--) {
	  if (u(k)) break;
	}
	if (k < 1) k = 1;
	for (i=k-1; i>=1; i--) {
	  rotate(r, qt, n, i, u(i), -u(i+1));
	  if (u(i) == 0.0) u(i) = fabs(u(i+1));
		else if (fabs(u(i)) > fabs(u(i+1)))
			u(i) = fabs(u(i)) * sqrt(1.0 + sqr(u(i+1) / u(i)));
		else u(i) = fabs(u(i+1)) * sqrt(1.0 + sqr(u(i) / u(i+1)));
	}
	for (j=1; j<=n; j++) r(1, j) += u(1) * v(j);
	for (i=1; i<k; i++)
	  rotate(r, qt, n, i, r(i, i), -r(i+1, i));
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

void rsolv(const DoubleMatrix & a, int n, const DoubleVector & d, 
	   DoubleVector & b) {
	int i,j;
	double sum;

	b(n) /= d(n);
	for (i=n-1; i>=1; i--) {
	  for (sum=0.0, j=i+1; j<=n; j++) sum += a(i, j) * b(j);
	  b(i)=(b(i) - sum) / d(i);
	}
}
