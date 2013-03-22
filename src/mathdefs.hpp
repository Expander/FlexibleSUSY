// -*- c++ -*-

#ifndef mathdefs_h
#define mathdefs_h


#include <cmath>
#include <complex>
#include <limits>
#include <iostream>		// required by <tvmet/Matrix.h>

#include <boost/numeric/ublas/vector.hpp>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>

typedef double Real;
typedef std::complex<Real> Comp;

typedef boost::numeric::ublas::vector<Real> BRVec;

typedef tvmet::Matrix<Real, 2, 2> RM22;
typedef tvmet::Matrix<Real, 3, 3> RM33;
typedef tvmet::Matrix<Real, 4, 4> RM44;
typedef tvmet::Matrix<Real, 6, 6> RM66;
typedef tvmet::Matrix<Comp, 2, 2> CM22;
typedef tvmet::Matrix<Comp, 3, 3> CM33;
typedef tvmet::Matrix<Comp, 4, 4> CM44;
typedef tvmet::Matrix<Comp, 6, 6> CM66;
typedef tvmet::Vector<Real, 2>	  RVe2;
typedef tvmet::Vector<Real, 3>	  RVe3;
typedef tvmet::Vector<Real, 4>	  RVe4;
typedef tvmet::Vector<Real, 6>	  RVe6;


const double pi = M_PI;
const double r2 = M_SQRT2;
const Real epsilon = std::numeric_limits<Real>::epsilon();


struct Interval {
    Real ini;			// initial value
    Real fin;			// final value
    size_t nsteps;		// number of steps
};


template<unsigned p, bool even_p, class T>
struct _power;

template<unsigned p, class T>
struct _power<p, true, T> {
    static T f(T x) {
	return _power<p/2, (p/2) % 2 == 0, T>::f(x) *
	       _power<p/2, (p/2) % 2 == 0, T>::f(x);
    }
};

template<unsigned p, class T>
struct _power<p, false, T> {
    static T f(T x) {
	return _power<p-1, (p-1) % 2 == 0, T>::f(x) * x;
    }
};

template<class T>
struct _power<1, false, T> {
    static T f(T x) {
	return x;
    }
};

template<class T>
struct _power<0, true, T> {
    static T f(T x) {
	return 1;
    }
};

template<unsigned p, class T> inline T Pow(T x)
{
    return _power<p, p % 2 == 0, T>::f(x);
}

template<class T> inline T Pow2(T x) { return Pow<2>(x); }
template<class T> inline T Pow3(T x) { return Pow<3>(x); }
template<class T> inline T Pow4(T x) { return Pow<4>(x); }
template<class T> inline T Pow5(T x) { return Pow<5>(x); }
template<class T> inline T sqr (T x) { return Pow<2>(x); }
template<class T> inline T cube(T x) { return Pow<3>(x); }

template<class T> inline T sign(T x)
{
    if (x == 0) return 0;
    else if (x > 0) return 1;
    else return -1;
}


#endif // mathdefs_h
