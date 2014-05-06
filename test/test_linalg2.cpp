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
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/placeholders.hpp>

using namespace std;
using namespace Eigen;
using namespace flexiblesusy;

template<class S_, int N_,
	 void fxn_(const Matrix<S_, N_, N_>&,
		   Array<double, N_, 1>&,
		   Matrix<S_, N_, N_>&,
		   Matrix<S_, N_, N_>&),
	 bool check_ascending_order_ = false>
struct Test_svd {
    typedef S_ S;
    enum { N = N_ };
    enum { check_ascending_order = check_ascending_order_ };
    void fxn(const Matrix<S_, N_, N_>& m,
	     Array<double, N_, 1>& s,
	     Matrix<S_, N_, N_>& u,
	     Matrix<S_, N_, N_>& vh)
    { fxn_(m, s, u, vh); }
};

typedef boost::mpl::list<
    // use Eigen::JacobiSVD
    Test_svd<complex<double>, 2, svd>,
    Test_svd<complex<double>, 3, svd>,
    Test_svd<double	    , 2, svd>,
    Test_svd<double	    , 3, svd>,

    Test_svd<complex<double>, 2, reorder_svd, true>,
    Test_svd<complex<double>, 3, reorder_svd, true>,
    Test_svd<double	    , 2, reorder_svd, true>,
    Test_svd<double	    , 3, reorder_svd, true>,

    // use ZGESVD of LAPACK
    Test_svd<complex<double>, 4, svd>,
    Test_svd<complex<double>, 6, svd>,

    Test_svd<complex<double>, 4, reorder_svd, true>,
    Test_svd<complex<double>, 6, reorder_svd, true>,

    // use DGESVD of LAPACK
    Test_svd<double	    , 4, svd>,
    Test_svd<double	    , 6, svd>,

    Test_svd<double	    , 4, reorder_svd, true>,
    Test_svd<double	    , 6, reorder_svd, true>
> svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_svd, T, svd_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<S, N, N> u, vh;

    T().fxn(m, s, u, vh);	// following LAPACK convention
    Matrix<S, N, N> diag = u.adjoint() * m * vh.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);

    if (T::check_ascending_order)
	for (size_t i = 0; i < N-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);
}

template<class S_, int N_,
	 void fxn_(const Matrix<S_, N_, N_>&,
		   Array<double, N_, 1>&, Matrix<complex<double>, N_, N_>&),
	 bool check_ascending_order_ = false>
struct Test_diagonalize_symmetric {
    typedef S_ S;
    enum { N = N_ };
    enum { check_ascending_order = check_ascending_order_ };
    void fxn(const Matrix<S_, N_, N_>& m,
	     Array<double, N_, 1>& s, Matrix<complex<double>, N_, N_>& u)
    { fxn_(m, s, u); }
};

typedef boost::mpl::list<
    // use Eigen::JacobiSVD
    Test_diagonalize_symmetric<complex<double>, 2, diagonalize_symmetric>,
    Test_diagonalize_symmetric<complex<double>, 3, diagonalize_symmetric>,

    Test_diagonalize_symmetric
	<complex<double>, 2, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<complex<double>, 3, reorder_diagonalize_symmetric, true>,

    // use Eigen::SelfAdjointEigenSolver
    Test_diagonalize_symmetric<double, 6, diagonalize_symmetric>,

    Test_diagonalize_symmetric<double, 6, reorder_diagonalize_symmetric, true>,

    // use ZGESVD of LAPACK
    Test_diagonalize_symmetric<complex<double>, 4, diagonalize_symmetric>,
    Test_diagonalize_symmetric<complex<double>, 6, diagonalize_symmetric>,

    Test_diagonalize_symmetric
	<complex<double>, 4, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<complex<double>, 6, reorder_diagonalize_symmetric, true>
> diagonalize_symmetric_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_diagonalize_symmetric, T, diagonalize_symmetric_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u;

    T().fxn(m, s, u);
    Matrix<complex<double>, N, N> diag = u.adjoint() * m * u.conjugate();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);

    if (T::check_ascending_order)
	for (size_t i = 0; i < N-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);
}

template<class S_, int N_,
	 void fxn_(const Matrix<S_, N_, N_>&,
		   Array<double, N_, 1>&,
		   Matrix<S_, N_, N_>&)>
struct Test_diagonalize_hermitian {
    typedef S_ S;
    enum { N = N_ };
    void fxn(const Matrix<S_, N_, N_>& m,
	     Array<double, N_, 1>& w,
	     Matrix<S_, N_, N_>& z)
    { fxn_(m, w, z); }
};

typedef boost::mpl::list<
    // use Eigen::SelfAdjointEigenSolver
    Test_diagonalize_hermitian<complex<double>, 6, hermitian_eigen>,
    Test_diagonalize_hermitian<double	      , 6, hermitian_eigen>,

    // use ZHEEV of LAPACK
    Test_diagonalize_hermitian<complex<double>, 6, hermitian_lapack>,
    // use DSYEV of LAPACK
    Test_diagonalize_hermitian<double	      , 6, hermitian_lapack>
> diagonalize_hermitian_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_diagonalize_hermitian, T, diagonalize_hermitian_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

    Matrix<S, N, N> m = Matrix<S, N, N>::Random();
    m = ((m + m.adjoint())/2).eval();
    Array<double, N, 1> w;
    Matrix<S, N, N> z;

    T().fxn(m, w, z);		// following LAPACK convention
    Matrix<S, N, N> diag = z.adjoint() * m * z;

    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 1e-12);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(w[i] <= w[i+1]);
}

template<class S_, int N_>
struct Test_fs {
    typedef S_ S;
    enum { N = N_ };
};

typedef boost::mpl::list<
    Test_fs<complex<double>, 2>,
    Test_fs<complex<double>, 3>,
    Test_fs<complex<double>, 4>,
    Test_fs<complex<double>, 6>,
    Test_fs<double	   , 2>,
    Test_fs<double	   , 3>,
    Test_fs<double	   , 4>,
    Test_fs<double	   , 6>
> fs_svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_fs_svd, T, fs_svd_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

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

typedef boost::mpl::list<
    boost::mpl::int_<2>,
    boost::mpl::int_<3>,
    boost::mpl::int_<4>,
    boost::mpl::int_<6>
> casting_fs_svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_casting_fs_svd, T, casting_fs_svd_tests)
{
    const size_t N = T::value;

    Matrix<double, N, N> m = Matrix<double, N, N>::Random();
    Array<double, N, 1> s;
    Matrix<complex<double>, N, N> u, v;

    fs_svd(m, s, u, v);		// following SARAH convention
    Matrix<complex<double>, N, N> diag = u.conjugate() * m * v.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (size_t i = 0; i < N; i++)
	for (size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-14);

    for (size_t i = 0; i < N-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);
}

typedef boost::mpl::list<
    // use Eigen::JacobiSVD
    Test_fs<complex<double>, 2>,
    Test_fs<complex<double>, 3>,

    // use ZGESVD of LAPACK
    Test_fs<complex<double>, 4>,
    Test_fs<complex<double>, 6>,

    // use Eigen::SelfAdjointEigenSolver
    Test_fs<double	   , 6>
> fs_diagonalize_symmetric_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_fs_diagonalize_symmetric, T, fs_diagonalize_symmetric_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

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

using namespace boost::mpl::placeholders;

typedef boost::mpl::fold<
    boost::mpl::range_c<int, 0, 50>,
    boost::mpl::list<>,
    boost::mpl::push_front<
	boost::mpl::push_front<_1, Test_fs<complex<double>, 6> >,
	Test_fs<double, 6> >
>::type fs_diagonalize_hermitian_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_fs_diagonalize_hermitian, T, fs_diagonalize_hermitian_tests)
{
    typedef typename T::S S;
    const size_t N = T::N;

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
