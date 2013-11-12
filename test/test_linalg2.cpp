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

#include <cmath>
#include <complex>
#include "linalg2.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_linalg2

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace Eigen;
using namespace flexiblesusy;

template<class S, int N, decltype(reorder_svd<S, N, N>) *fxn>
void test_svd(bool check_ascending_order = false)
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<S, N, N> u, vh;

    fxn(m, s, u, vh);		// following LAPACK convention
    Matrix<S, N, N> diag = u.adjoint() * m * vh.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);

    if (check_ascending_order)
	for (size_t i = 0; i < N-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);
}

BOOST_AUTO_TEST_CASE(test_svd_eigen)
{
    // uses Eigen::JacobiSVD
    test_svd<complex<double>, 2, svd>();
    test_svd<complex<double>, 3, svd>();
    test_svd<double	    , 2, svd>();
    test_svd<double	    , 3, svd>();

    test_svd<complex<double>, 2, reorder_svd>(true);
    test_svd<complex<double>, 3, reorder_svd>(true);
    test_svd<double	    , 2, reorder_svd>(true);
    test_svd<double	    , 3, reorder_svd>(true);
}

BOOST_AUTO_TEST_CASE(test_svd_lapack)
{
    // uses ZGESVD of LAPACK
    test_svd<complex<double>, 4, svd>();
    test_svd<complex<double>, 6, svd>();

    test_svd<complex<double>, 4, reorder_svd>(true);
    test_svd<complex<double>, 6, reorder_svd>(true);

    // uses DGESVD of LAPACK
    test_svd<double	    , 4, svd>();
    test_svd<double	    , 6, svd>();

    test_svd<double	    , 4, reorder_svd>(true);
    test_svd<double	    , 6, reorder_svd>(true);
}

template<class S, int N,
	 void fxn(const Matrix<S, N, N>&,
		  Array<double, N, 1>&, Matrix<complex<double>, N, N>&)>
void test_diagonalize_symmetric(bool check_ascending_order = false)
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u;

    fxn(m, s, u);
    Matrix<complex<double>, N, N> diag = u.adjoint() * m * u.conjugate();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);

    if (check_ascending_order)
	for (size_t i = 0; i < N-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);
}

BOOST_AUTO_TEST_CASE(test_diagonalize_symmetric_eigen)
{
    // uses Eigen::JacobiSVD
    test_diagonalize_symmetric<complex<double>, 2, diagonalize_symmetric>();
    test_diagonalize_symmetric<complex<double>, 3, diagonalize_symmetric>();

    test_diagonalize_symmetric
	<complex<double>, 2, reorder_diagonalize_symmetric>(true);
    test_diagonalize_symmetric
	<complex<double>, 3, reorder_diagonalize_symmetric>(true);

    // uses Eigen::SelfAdjointEigenSolver
    test_diagonalize_symmetric<double, 6, diagonalize_symmetric>();

    test_diagonalize_symmetric<double, 6, reorder_diagonalize_symmetric>(true);
}

BOOST_AUTO_TEST_CASE(test_diagonalize_symmetric_lapack)
{
    // uses ZGESVD of LAPACK
    test_diagonalize_symmetric<complex<double>, 4, diagonalize_symmetric>();
    test_diagonalize_symmetric<complex<double>, 6, diagonalize_symmetric>();

    test_diagonalize_symmetric
	<complex<double>, 4, reorder_diagonalize_symmetric>(true);
    test_diagonalize_symmetric
	<complex<double>, 6, reorder_diagonalize_symmetric>(true);
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

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(w[i] <= w[i+1]);
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

template<class S, int N>
void check_fs_svd()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<S, N, N> u, v;

    fs_svd(m, s, u, v);		// following SARAH convention
    Matrix<S, N, N> diag = u.conjugate() * m * v.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);
}

BOOST_AUTO_TEST_CASE(test_fs_svd)
{
    check_fs_svd<complex<double>, 2>();
    check_fs_svd<complex<double>, 3>();
    check_fs_svd<complex<double>, 4>();
    check_fs_svd<complex<double>, 6>();
    check_fs_svd<double	   	, 2>();
    check_fs_svd<double	   	, 3>();
    check_fs_svd<double	   	, 4>();
    check_fs_svd<double	   	, 6>();
}

template<int N>
void check_casting_fs_svd()
{
    Matrix<double, N, N> m = Matrix<double, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u, v;

    fs_svd(m, s, u, v);	// following SARAH convention
    Matrix<complex<double>, N, N> diag = u.conjugate() * m * v.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);
}

BOOST_AUTO_TEST_CASE(test_casting_fs_svd)
{
    check_casting_fs_svd<2>();
    check_casting_fs_svd<3>();
    check_casting_fs_svd<4>();
    check_casting_fs_svd<6>();
}

template<class S, int N>
void check_fs_diagonalize_symmetric()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u;

    fs_diagonalize_symmetric(m, s, u);
    Matrix<complex<double>, N, N> diag = u.conjugate() * m * u.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);
}

BOOST_AUTO_TEST_CASE(test_fs_diagonalize_symmetric)
{
    // uses Eigen::JacobiSVD
    check_fs_diagonalize_symmetric<complex<double>, 2>();
    check_fs_diagonalize_symmetric<complex<double>, 3>();

    // uses ZGESVD of LAPACK
    check_fs_diagonalize_symmetric<complex<double>, 4>();
    check_fs_diagonalize_symmetric<complex<double>, 6>();

    // uses Eigen::SelfAdjointEigenSolver
    check_fs_diagonalize_symmetric<double	  , 6>();
}

template<class S, int N>
void check_fs_diagonalize_hermitian()
{
    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.adjoint())/2).eval();
    Array<double, N, 1> w;
    Matrix<S, N, N> z;

    fs_diagonalize_hermitian(m, w, z); // following SARAH convention
    Matrix<S, N, N> diag = z * m * z.adjoint();

    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 1e-11);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(abs(w[i]) <= abs(w[i+1]));
}

BOOST_AUTO_TEST_CASE(test_fs_diagonalize_hermitian)
{
    for (size_t n = 50; n; n--) {
	check_fs_diagonalize_hermitian<complex<double>, 6>();
	check_fs_diagonalize_hermitian<double	      , 6>();
    }
}
