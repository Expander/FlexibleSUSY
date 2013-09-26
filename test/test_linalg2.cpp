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

#include <complex>
#include "linalg2.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_linalg2

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace Eigen;
using namespace flexiblesusy;

template<class S, int N>
void test_svd()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<S, N, N> u, vh;

    svd(m, s, u, vh);		// following LAPACK convention
    Matrix<S, N, N> diag = u.adjoint() * m * vh.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);
}

BOOST_AUTO_TEST_CASE(test_svd_eigen)
{
    // uses Eigen::JacobiSVD
    test_svd<complex<double>, 2>();
    test_svd<complex<double>, 3>();
    test_svd<double	    , 2>();
    test_svd<double	    , 3>();
}

BOOST_AUTO_TEST_CASE(test_svd_lapack)
{
    // uses ZGESVD of LAPACK
    test_svd<complex<double>, 4>();
    test_svd<complex<double>, 6>();
    // uses DGESVD of LAPACK
    test_svd<double	    , 4>();
    test_svd<double	    , 6>();
}

template<class S, int N>
void test_diagonalize_symmetric()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u;

    diagonalize_symmetric(m, s, u);
    Matrix<complex<double>, N, N> diag = u.adjoint() * m * u.conjugate();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-13);
}

BOOST_AUTO_TEST_CASE(test_diagonalize_symmetric_eigen)
{
    // uses Eigen::JacobiSVD
    test_diagonalize_symmetric<complex<double>, 2>();
    test_diagonalize_symmetric<complex<double>, 3>();
    // uses Eigen::SelfAdjointEigenSolver
    test_diagonalize_symmetric<double	      , 6>();
}

BOOST_AUTO_TEST_CASE(test_diagonalize_symmetric_lapack)
{
    // uses ZGESVD of LAPACK
    test_diagonalize_symmetric<complex<double>, 4>();
    test_diagonalize_symmetric<complex<double>, 6>();
}

template<class S, int N, decltype(diagonalize_hermitian<S, N>) *fxn>
void test_diagonalize_hermitian()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.adjoint())/2).eval();
    Array<double, N, 1> w;
    Matrix<S, N, N> z;

    fxn(m, w, z);		// following LAPACK convention
    Matrix<S, N, N> diag = z.adjoint() * m * z;

    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 1e-12);
}

BOOST_AUTO_TEST_CASE(test_diagonalize_hermitian_eigen)
{
    // uses Eigen::SelfAdjointEigenSolver
    test_diagonalize_hermitian<complex<double>, 6, hermitian_eigen>();
    test_diagonalize_hermitian<double	      , 6, hermitian_eigen>();
}

BOOST_AUTO_TEST_CASE(test_diagonalize_hermitian_lapack)
{
    // uses ZHEEV of LAPACK
    test_diagonalize_hermitian<complex<double>, 6, hermitian_lapack>();
    // uses DSYEV of LAPACK
    test_diagonalize_hermitian<double	      , 6, hermitian_lapack>();
}
