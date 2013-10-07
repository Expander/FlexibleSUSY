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

using namespace flexiblesusy;

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
