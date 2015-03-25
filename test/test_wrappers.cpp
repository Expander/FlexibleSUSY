// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_wrappers

#include <random>
#include <complex>
#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"
#include "diagonalization.hpp"
#include "stopwatch.hpp"
#include <boost/lexical_cast.hpp>

#if defined(__CYGWIN__) || defined(__FreeBSD__)
// std::to_string is missing on Cygwin and FreeBSD
// see http://stackoverflow.com/questions/22571838/gcc-4-8-1-stdto-string-error
#   undef  ENABLE_STD_TO_STRING
#else
#   define ENABLE_STD_TO_STRING 1
#endif

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_Delta )
{
   BOOST_CHECK_EQUAL(Delta(0,0), 1);
   BOOST_CHECK_EQUAL(Delta(1,1), 1);
   BOOST_CHECK_EQUAL(Delta(2,2), 1);
   BOOST_CHECK_EQUAL(Delta(3,3), 1);

   BOOST_CHECK_EQUAL(Delta(-1,-1), 1);
   BOOST_CHECK_EQUAL(Delta(-2,-2), 1);
   BOOST_CHECK_EQUAL(Delta(-3,-3), 1);

   BOOST_CHECK_EQUAL(Delta(0,1), 0);
   BOOST_CHECK_EQUAL(Delta(1,0), 0);
   BOOST_CHECK_EQUAL(Delta(0,2), 0);
   BOOST_CHECK_EQUAL(Delta(2,0), 0);

   BOOST_CHECK_EQUAL(Delta(0,-1), 0);
   BOOST_CHECK_EQUAL(Delta(-1,0), 0);
   BOOST_CHECK_EQUAL(Delta(0,-2), 0);
   BOOST_CHECK_EQUAL(Delta(-2,0), 0);
}

using namespace std;

DoubleMatrix random_real_matrix(int n, int m)
{
    static default_random_engine generator;
    static uniform_real_distribution<> o1(-3, 3);

    DoubleMatrix r(n, m);
    for (int i = 1; i <= n; i++)
	for (int j = 1; j <= n; j++)
	    r(i, j) = o1(generator);
    return r;
}

BOOST_AUTO_TEST_CASE(test_svd)
{
    for (int n = 2; n <= 6; n++) {
	DoubleMatrix  m(n,n);
	ComplexMatrix u(n,n);
	ComplexMatrix v(n,n);
	DoubleVector  s(n);
	ComplexMatrix diag(n,n);

	for (int count = 100; count; count--) {
	    m = random_real_matrix(n,n);
	    if (n == 2)
		Diagonalize2by2(m, u, v, s);
	    else
		Diagonalize(m, u, v, s);
	    diag = u.complexConjugate() * m * v.hermitianConjugate();

	    for (int i = 1; i <= s.displayEnd(); i++)
		BOOST_CHECK(s(i) >= 0);
	    for (int i = 1; i <= diag.displayCols(); i++)
		for (int j = 1; j <= diag.displayRows(); j++)
		    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j?s(i):0)), 1e-13);
	}
    }
}

BOOST_AUTO_TEST_CASE(test_symmetric)
{
    for (int n = 2; n <= 6; n++) {
	DoubleMatrix  m(n,n);
	ComplexMatrix u(n,n);
	DoubleVector  s(n);
	ComplexMatrix diag(n,n);

	for (int count = 100; count; count--) {
	    m = random_real_matrix(n,n);
	    m.symmetrise();
	    if (n == 2)
		Diagonalize2by2(m, u, s);
	    else
		Diagonalize(m, u, s);
	    diag = u.complexConjugate() * m * u.hermitianConjugate();

	    for (int i = 1; i <= s.displayEnd(); i++)
		BOOST_CHECK(s(i) >= 0);
	    for (int i = 1; i <= diag.displayCols(); i++)
		for (int j = 1; j <= diag.displayRows(); j++)
		    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j?s(i):0)), 1e-13);
	}
    }
}

BOOST_AUTO_TEST_CASE(test_Diag)
{
   Eigen::Matrix<double,3,3> m;
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         m(i,k) = (i+1) * (k+1);

   Eigen::Matrix<double,3,3> diag(Diag(m));

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++) {
         if (i == k)
            BOOST_CHECK_PREDICATE(std::not_equal_to<double>(), (diag(i,k))(0.));
         else
            BOOST_CHECK_EQUAL(diag(i,k), 0.);
      }
}

template <typename T>
std::string ToString_sstream(T a)
{
   std::ostringstream ostr;
   ostr << a;
   return ostr.str();
}

#ifdef ENABLE_STD_TO_STRING
template <typename T>
std::string ToString_to_string(T a)
{
   return std::to_string(a);
}
#endif

template <typename T>
std::string ToString_sprintf(T a)
{
   static const unsigned buf_length = 20;
   char buf[buf_length];
   snprintf(buf, buf_length, "%i");
   return std::string(buf);
}

template <typename T>
std::string ToString_lexical_cast(T a)
{
   return boost::lexical_cast<std::string>(a);
}

#define MEASURE(type,number,iterations)                            \
   do {                                                            \
      Stopwatch stopwatch;                                         \
      double time = 0.;                                            \
      for (int i = 0; i < iterations; i++) {                       \
         stopwatch.start();                                        \
         ToString_##type(number);                                  \
         stopwatch.stop();                                         \
         time += stopwatch.get_time_in_seconds();                  \
      }                                                            \
      BOOST_MESSAGE("ToString via " #type ": " << time << " s");   \
   } while (0)

BOOST_AUTO_TEST_CASE(test_ToString)
{
   const int number = 123456;
   const int number_of_iterations = 1000000;

   MEASURE(sstream     , number, number_of_iterations);
#ifdef ENABLE_STD_TO_STRING
   MEASURE(to_string   , number, number_of_iterations);
#endif
   MEASURE(sprintf     , number, number_of_iterations);
   MEASURE(lexical_cast, number, number_of_iterations);
}

BOOST_AUTO_TEST_CASE(test_calculate_singlet_mass)
{
   double mass;
   std::complex<double> phase;

   mass = calculate_singlet_mass(100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_singlet_mass(-100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_SMALL(phase.real(), 1e-15);
   BOOST_CHECK_EQUAL(phase.imag(), 1.);

   mass = calculate_singlet_mass(std::complex<double>(100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_singlet_mass(std::complex<double>(-100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_SMALL(phase.real(), 1e-15);
   BOOST_CHECK_EQUAL(phase.imag(), 1.);

   mass = calculate_singlet_mass(std::complex<double>(1.,1.), phase);
   BOOST_CHECK_EQUAL(mass, sqrt(2.));
   BOOST_CHECK_EQUAL(std::abs(phase), 1.);
   BOOST_CHECK_CLOSE(std::arg(phase), 0.5 * Pi/4., 1e-13);
}
