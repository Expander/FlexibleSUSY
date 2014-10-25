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

#ifndef linalg2_hpp
#define linalg2_hpp

#include <limits>
#include <cmath>
#include <complex>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

namespace flexiblesusy {

#define MAX_(i, j) (((i) > (j)) ? (i) : (j))
#define MIN_(i, j) (((i) < (j)) ? (i) : (j))

template<class Scalar, int M, int N>
void svd_eigen
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v)
{
    Eigen::JacobiSVD<Eigen::Matrix<Scalar, M, N> >
	svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
    s = svd.singularValues();
    u = svd.matrixU();
    v = svd.matrixV();
}

template<class Scalar, int M, int N>
void svd_eigen
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u)
{
    Eigen::JacobiSVD<Eigen::Matrix<Scalar, M, N> > svd(m, Eigen::ComputeFullU);
    s = svd.singularValues();
    u = svd.matrixU();
}

template<class Scalar, int N>
void hermitian_eigen
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar,N,N> > es(m);
    w = es.eigenvalues();
    z = es.eigenvectors();
}

extern "C" void zgesvd_
(const char& JOBU, const char& JOBVT, const int& M, const int& N,
 std::complex<double> *A, const int& LDA, double *S, std::complex<double> *U,
 const int& LDU, std::complex<double> *VT, const int& LDVT,
 std::complex<double> *WORK, const int& LWORK, double *RWORK, int& INFO);

extern "C" void dgesvd_
(const char& JOBU, const char& JOBVT, const int& M, const int& N,
 double *A, const int& LDA, double *S, double *U,
 const int& LDU, double *VT, const int& LDVT,
 double *WORK, const int& LWORK, int& INFO);

extern "C" void zheev_
(const char& JOBZ, const char& UPLO, const int& N, std::complex<double> *A,
 const int& LDA, double *W, std::complex<double> *WORK, const int& LWORK,
 double *RWORK, int& INFO);

extern "C" void dsyev_
(const char& JOBZ, const char& UPLO, const int& N, double *A,
 const int& LDA, double *W, double *WORK, const int& LWORK,
 int& INFO);

extern "C" void ddisna_
(const char& JOB, const int& M, const int& N, const double *D, double *SEP,
 int& INFO);

#define def_svd_lapack(t, f, ...)					\
template<int M, int N>							\
void svd_lapack								\
(const Eigen::Matrix<t, M, N>& m,					\
 Eigen::Array<double, MIN_(M, N), 1>& s,				\
 Eigen::Matrix<t, M, M>& u,						\
 Eigen::Matrix<t, N, N>& vh)						\
{									\
    call_lapack_svd(t, f, 'A', vh.data(), __VA_ARGS__);			\
}									\
									\
template<int M, int N>							\
void svd_lapack								\
(const Eigen::Matrix<t, M, N>& m,					\
 Eigen::Array<double, MIN_(M, N), 1>& s,				\
 Eigen::Matrix<t, M, M>& u)						\
{									\
    call_lapack_svd(t, f, 'N', 0, __VA_ARGS__);				\
}

#define call_lapack_svd(t, f, jobvt, vt, ...)				\
    const     char JOBU  = 'A';						\
    const     char JOBVT = (jobvt);					\
    Eigen::Matrix<t, M, N> A = m;					\
    const     int LDA   = M;						\
    const     int LDU   = M;						\
    const     int LDVT  = N;						\
    const     int LWORK = get_lwork(__VA_ARGS__,);			\
    Eigen::Array<t, LWORK, 1> WORK;					\
    decl_rwork(__VA_ARGS__);						\
    int INFO;								\
    f(JOBU, JOBVT, M, N, A.data(), LDA, s.data(), u.data(), LDU,	\
      (vt), LDVT, WORK.data(), LWORK, put_rwork(__VA_ARGS__) INFO);

#define def_hermitian_lapack(s, f, ...)					\
template<int N>								\
void hermitian_lapack							\
(const Eigen::Matrix<s, N, N>& m,					\
 Eigen::Array<double, N, 1>& w,						\
 Eigen::Matrix<s, N, N>& z)						\
{									\
    const     char JOBZ = 'V';						\
    const     char UPLO = 'L';						\
    Eigen::Matrix<s, N, N> A = m;					\
    const     int LDA   = N;						\
    const     int LWORK = get_lwork(__VA_ARGS__,);			\
    Eigen::Array<s, LWORK, 1> WORK;					\
    decl_rwork(__VA_ARGS__);						\
    int INFO;								\
    f(JOBZ, UPLO, N, A.data(), LDA, w.data(), WORK.data(), LWORK,	\
      put_rwork(__VA_ARGS__) INFO);					\
    z = A;								\
}

#define get_lwork(lwork, ...) (lwork)

#define get_rwork_macro(_1, _2, name, ...) name

#define nop_(_1)

#define do_decl_rwork(_1, lrwork) Eigen::Array<double, (lrwork), 1> RWORK

#define decl_rwork(...) \
    get_rwork_macro(__VA_ARGS__, do_decl_rwork, nop_,)(__VA_ARGS__)

#define do_put_rwork(_1, _2) RWORK.data(),

#define put_rwork(...) \
    get_rwork_macro(__VA_ARGS__, do_put_rwork, nop_,)(__VA_ARGS__)

def_svd_lapack(std::complex<double>, zgesvd_, 3*MAX_(M,N), 5*MIN_(M,N))
def_svd_lapack(double, dgesvd_, MAX_(3*MIN_(M,N)+MAX_(M,N),5*MIN_(M,N)))

def_hermitian_lapack(std::complex<double>, zheev_, 2*N-1, 3*N-2)
def_hermitian_lapack(double, dsyev_, 3*N-1)


// ZGESVD of ATLAS seems to be faster than Eigen::JacobiSVD for M, N >= 4

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u * sigma * vh    // LAPACK convention
 *
 * and `(s >= 0).all()`.  Elements of s are in descending order.  The
 * above decomposition can be put in the form
 *
 *     m == u * s.matrix().asDiagonal() * vh
 *
 * if `M == N`.
 *
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] vh     N-by-N unitary matrix
 */
template<class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh)
{
    svd_lapack(m, s, u, vh);
}

template<class Scalar>
void svd
(const Eigen::Matrix<Scalar, 3, 3>& m,
 Eigen::Array<double, 3, 1>& s,
 Eigen::Matrix<Scalar, 3, 3>& u,
 Eigen::Matrix<Scalar, 3, 3>& vh)
{
    svd_eigen(m, s, u, vh);
    vh.adjointInPlace();
}

template<class Scalar>
void svd
(const Eigen::Matrix<Scalar, 2, 2>& m,
 Eigen::Array<double, 2, 1>& s,
 Eigen::Matrix<Scalar, 2, 2>& u,
 Eigen::Matrix<Scalar, 2, 2>& vh)
{
    svd_eigen(m, s, u, vh);
    vh.adjointInPlace();
}

template<class Scalar, int M, int N>
void svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double *s_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *v_errbd = 0)
{
    svd(m, s, u, vh);

    // see http://www.netlib.org/lapack/lug/node96.html
    if (!s_errbd) return;
    const double EPSMCH = std::numeric_limits<double>::epsilon();
    *s_errbd = EPSMCH * s[0];

    Eigen::Array<double, MIN_(M, N), 1> RCOND;
    int INFO;
    if (u_errbd) {
	ddisna_('L', M, N, s.data(), RCOND.data(), INFO);
	u_errbd->fill(*s_errbd);
	*u_errbd /= RCOND;
    }
    if (v_errbd) {
	ddisna_('R', M, N, s.data(), RCOND.data(), INFO);
	v_errbd->fill(*s_errbd);
	*v_errbd /= RCOND;
    }
}

template<class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u)
{
    svd_lapack(m, s, u);
}

template<class Scalar>
void svd
(const Eigen::Matrix<Scalar, 3, 3>& m,
 Eigen::Array<double, 3, 1>& s,
 Eigen::Matrix<Scalar, 3, 3>& u)
{
    svd_eigen(m, s, u);
}

template<class Scalar>
void svd
(const Eigen::Matrix<Scalar, 2, 2>& m,
 Eigen::Array<double, 2, 1>& s,
 Eigen::Matrix<Scalar, 2, 2>& u)
{
    svd_eigen(m, s, u);
}

template<class Scalar, int M, int N>
void svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 double *s_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *u_errbd = 0)
{
    svd(m, s, u);

    // see http://www.netlib.org/lapack/lug/node96.html
    if (!s_errbd) return;
    const double EPSMCH = std::numeric_limits<double>::epsilon();
    *s_errbd = EPSMCH * s[0];

    Eigen::Array<double, MIN_(M, N), 1> RCOND;
    int INFO;
    if (u_errbd) {
	ddisna_('L', M, N, s.data(), RCOND.data(), INFO);
	u_errbd->fill(*s_errbd);
	*u_errbd /= RCOND;
    }
}

/**
 * Same as svd(m, s, u, vh) except that an approximate error bound for
 * the singular values is returned.  The error bound is estimated
 * following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of svd(m, s, u, vh) for the other parameters.
 */
template<class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double& s_errbd)
{
    svd_errbd(m, s, u, vh, &s_errbd);
}

/**
 * Same as svd(m, s, u, vh, s_errbd) except that approximate error
 * bounds for the singular vectors are returned.  The error bounds are
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of svd(m, s, u, vh, s_errbd) for the other
 * parameters.
 */
template<class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double& s_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& v_errbd)
{
    svd_errbd(m, s, u, vh, &s_errbd, &u_errbd, &v_errbd);
}

// Eigen::SelfAdjointEigenSolver seems to be faster than ZHEEV of ATLAS

/**
 * Diagonalizes N-by-N hermitian matrix m so that
 *
 *     m == z * w.matrix().asDiagonal() * z.adjoint()
 *
 * Elements of w are in ascending order.
 *
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 * @param[out] z      N-by-N unitary matrix
 */
template<class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z)
{
    hermitian_eigen(m, w, z);
}

template<class Scalar, int N>
void diagonalize_hermitian_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double *w_errbd = 0,
 Eigen::Array<double, N, 1> *z_errbd = 0)
{
    diagonalize_hermitian(m, w, z);

    // see http://www.netlib.org/lapack/lug/node89.html
    if (!w_errbd) return;
    const double EPSMCH = std::numeric_limits<double>::epsilon();
    double mnorm = std::max(std::abs(w[0]), std::abs(w[N-1]));
    *w_errbd = EPSMCH * mnorm;

    if (!z_errbd) return;
    Eigen::Array<double, N, 1> RCONDZ;
    int INFO;
    ddisna_('E', N, N, w.data(), RCONDZ.data(), INFO);
    z_errbd->fill(*w_errbd);
    *z_errbd /= RCONDZ;
}

/**
 * Same as diagonalize_hermitian(m, w, z) except that an approximate
 * error bound for the eigenvalues is returned.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of diagonalize_hermitian(m, w, z) for the
 * other parameters.
 */
template<class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double& w_errbd)
{
    diagonalize_hermitian_errbd(m, w, z, &w_errbd);
}

/**
 * Same as diagonalize_hermitian(m, w, z, w_errbd) except that
 * approximate error bounds for the eigenvectors are returned.  The
 * error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] z_errbd array of approximate error bounds for z
 *
 * See the documentation of diagonalize_hermitian(m, w, z, w_errbd)
 * for the other parameters.
 */
template<class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double& w_errbd,
 Eigen::Array<double, N, 1>& z_errbd)
{
    diagonalize_hermitian_errbd(m, w, z, &w_errbd, &z_errbd);
}

struct RephaseOp {
    std::complex<double> operator() (const std::complex<double>& z) const
	{ return std::polar(1.0, std::arg(z)/2); }
};

template<int N>
void diagonalize_symmetric_errbd
(const Eigen::Matrix<std::complex<double>, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double *s_errbd = 0,
 Eigen::Array<double, N, 1> *u_errbd = 0)
{
    svd_errbd(m, s, u, s_errbd, u_errbd);
    Eigen::Array<std::complex<double>, N, 1> diag =
	(u.adjoint() * m * u.conjugate()).diagonal();
    u *= diag.unaryExpr(RephaseOp()).matrix().asDiagonal();
}

/**
 * Diagonalizes N-by-N complex symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Elements of s are in descending order.
 *
 * @tparam     N number of rows and columns in m and u
 * @param[in]  m N-by-N complex symmetric matrix to be decomposed
 * @param[out] s array of length N to contain singular values
 * @param[out] u N-by-N complex unitary matrix
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<double>, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u)
{
    diagonalize_symmetric_errbd(m, s, u);
}

/**
 * Same as diagonalize_symmetric(m, s, u) except that an approximate
 * error bound for the singular values is returned.  The error bound
 * is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<double>, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, u, &s_errbd);
}

/**
 * Same as diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<double>, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd,
 Eigen::Array<double, N, 1>& u_errbd)
{
    diagonalize_symmetric_errbd(m, s, u, &s_errbd, &u_errbd);
}

struct FlipSignOp {
    std::complex<double> operator() (const std::complex<double>& z) const {
	return z.real() < 0 ? std::complex<double>(0.0,1.0) :
	    std::complex<double>(1.0,0.0);
    }
};

template<int N>
void diagonalize_symmetric_errbd
(const Eigen::Matrix<double, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double *s_errbd = 0,
 Eigen::Array<double, N, 1> *u_errbd = 0)
{
    Eigen::Matrix<double, N, N> z;
    diagonalize_hermitian_errbd(m, s, z, s_errbd, u_errbd);
    // see http://forum.kde.org/viewtopic.php?f=74&t=62606
    u = z * s.template cast<std::complex<double> >().
	unaryExpr(FlipSignOp()).matrix().asDiagonal();
    s = s.abs();
}

/**
 * Diagonalizes N-by-N real symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Order of elements of s is *unspecified*.
 *
 * @tparam     N number of rows and columns of m
 * @param[in]  m N-by-N real symmetric matrix to be decomposed
 * @param[out] s array of length N to contain singular values
 * @param[out] u N-by-N complex unitary matrix
 *
 * @note Use diagonalize_hermitian() unless sign of `s[i]` matters.
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<double, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u)
{
    diagonalize_symmetric_errbd(m, s, u);
}

/**
 * Same as diagonalize_symmetric(m, s, u) except that an approximate
 * error bound for the singular values is returned.  The error bound
 * is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<double, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, u, &s_errbd);
}

/**
 * Same as diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<int N>
void diagonalize_symmetric
(const Eigen::Matrix<double, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd,
 Eigen::Array<double, N, 1>& u_errbd)
{
    diagonalize_symmetric_errbd(m, s, u, &s_errbd, &u_errbd);
}

template<class Scalar, int M, int N>
void reorder_svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double *s_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *v_errbd = 0)
{
    svd_errbd(m, s, u, vh, s_errbd, u_errbd, v_errbd);
    s.reverseInPlace();
    Eigen::PermutationMatrix<MIN_(M, N)> p;
    p.setIdentity();
    p.indices().reverseInPlace();
    u              *= p;
    vh.transpose() *= p;
    if (u_errbd) u_errbd->reverseInPlace();
    if (v_errbd) v_errbd->reverseInPlace();
}

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u * sigma * vh    // LAPACK convention
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u * s.matrix().asDiagonal() * vh
 *
 * if `M == N`.
 *
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] vh     N-by-N unitary matrix
 */
template<class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh)
{
    reorder_svd_errbd(m, s, u, vh);
}

/**
 * Same as reorder_svd(m, s, u, vh) except that an approximate error
 * bound for the singular values is returned.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_svd(m, s, u, vh) for the other
 * parameters.
 */
template<class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double& s_errbd)
{
    reorder_svd_errbd(m, s, u, vh, &s_errbd);
}

/**
 * Same as reorder_svd(m, s, u, vh, s_errbd) except that approximate
 * error bounds for the singular vectors are returned.  The error
 * bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of reorder_svd(m, s, u, vh, s_errbd) for the
 * other parameters.
 */
template<class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 double& s_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& v_errbd)
{
    reorder_svd_errbd(m, s, u, vh, &s_errbd, &u_errbd, &v_errbd);
}

template<int N>
void reorder_diagonalize_symmetric_errbd
(const Eigen::Matrix<std::complex<double>, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double *s_errbd = 0,
 Eigen::Array<double, N, 1> *u_errbd = 0)
{
    diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    s.reverseInPlace();
    u = u.rowwise().reverse().eval();
    if (u_errbd) u_errbd->reverseInPlace();
}

template<int N>
struct Compare {
    Compare(const Eigen::Array<double, N, 1>& s_) : s(s_) {}
    bool operator() (int i, int j) { return s[i] < s[j]; }
    const Eigen::Array<double, N, 1>& s;
};

template<int N>
void reorder_diagonalize_symmetric_errbd
(const Eigen::Matrix<double, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double *s_errbd = 0,
 Eigen::Array<double, N, 1> *u_errbd = 0)
{
    diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    Eigen::PermutationMatrix<N> p;
    p.setIdentity();
    std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
	      Compare<N>(s));
#if EIGEN_VERSION_AT_LEAST(3,1,4)
    s.matrix().transpose() *= p;
    if (u_errbd) u_errbd->matrix().transpose() *= p;
#else
    Eigen::Map<Eigen::Matrix<double, N, 1> >(s.data()).transpose() *= p;
    if (u_errbd)
	Eigen::Map<Eigen::Matrix<double, N, 1> >(u_errbd->data()).transpose()
	    *= p;
#endif
    u *= p;
}

/**
 * Diagonalizes N-by-N symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 * @param[out] u      N-by-N complex unitary matrix
 */
template<class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u)
{
    reorder_diagonalize_symmetric_errbd(m, s, u);
}

/**
 * Same as reorder_diagonalize_symmetric(m, s, u) except that an
 * approximate error bound for the singular values is returned.  The
 * error bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_diagonalize_symmetric(m, s, u) for
 * the other parameters.
 */
template<class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd)
{
    reorder_diagonalize_symmetric_errbd(m, s, u, &s_errbd);
}

/**
 * Same as reorder_diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of reorder_diagonalize_symmetric(m, s, u,
 * s_errbd) for the other parameters.
 */
template<class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd,
 Eigen::Array<double, N, 1>& u_errbd)
{
    reorder_diagonalize_symmetric_errbd(m, s, u, &s_errbd, &u_errbd);
}

template<class Scalar, int M, int N>
void fs_svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v,
 double *s_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<double, MIN_(M, N), 1> *v_errbd = 0)
{
    reorder_svd_errbd(m, s, u, v, s_errbd, u_errbd, v_errbd);
    u.transposeInPlace();
}

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u.transpose() * sigma * v
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * v
 *
 * if `M == N`.
 *
 * @tparam     Scalar type of elements of m, u, and v
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] v      N-by-N unitary matrix
 */
template<class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v)
{
    fs_svd_errbd(m, s, u, v);
}

/**
 * Same as fs_svd(m, s, u, v) except that an approximate error bound
 * for the singular values is returned.  The error bound is estimated
 * following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_svd(m, s, u, v) for the other
 * parameters.
 */
template<class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v,
 double& s_errbd)
{
    fs_svd_errbd(m, s, u, v, &s_errbd);
}

/**
 * Same as fs_svd(m, s, u, v, s_errbd) except that approximate error
 * bounds for the singular vectors are returned.  The error bounds are
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of fs_svd(m, s, u, v, s_errbd) for the other
 * parameters.
 */
template<class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v,
 double& s_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& v_errbd)
{
    fs_svd_errbd(m, s, u, v, &s_errbd, &u_errbd, &v_errbd);
}

/**
 * Singular value decomposition of M-by-N *real* matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u.transpose() * sigma * v
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * v
 *
 * if `M == N`.
 *
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N *real* matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M *complex* unitary matrix
 * @param[out] v      N-by-N *complex* unitary matrix
 *
 * @note This is a convenience overload for the case where the type of
 * u and v (complex) differs from that of m (real).  Mathematically,
 * real u and v are enough to accommodate SVD of any real m.
 */
template<int M, int N>
void fs_svd
(const Eigen::Matrix<double, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<double>, M, M>& u,
 Eigen::Matrix<std::complex<double>, N, N>& v)
{
    fs_svd(m.template cast<std::complex<double> >().eval(), s, u, v);
}

/**
 * Same as fs_svd(m, s, u, v) except that an approximate error bound
 * for the singular values is returned.  The error bound is estimated
 * following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_svd(m, s, u, v) for the other
 * parameters.
 */
template<int M, int N>
void fs_svd
(const Eigen::Matrix<double, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<double>, M, M>& u,
 Eigen::Matrix<std::complex<double>, N, N>& v,
 double& s_errbd)
{
    fs_svd(m.template cast<std::complex<double> >().eval(), s, u, v, s_errbd);
}

/**
 * Same as fs_svd(m, s, u, v, s_errbd) except that approximate error
 * bounds for the singular vectors are returned.  The error bounds are
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of fs_svd(m, s, u, v, s_errbd) for the other
 * parameters.
 */
template<int M, int N>
void fs_svd
(const Eigen::Matrix<double, M, N>& m,
 Eigen::Array<double, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<double>, M, M>& u,
 Eigen::Matrix<std::complex<double>, N, N>& v,
 double& s_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<double, MIN_(M, N), 1>& v_errbd)
{
    fs_svd(m.template cast<std::complex<double> >().eval(), s, u, v,
	   s_errbd, u_errbd, v_errbd);
}

template<class Scalar, int N>
void fs_diagonalize_symmetric_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double *s_errbd = 0,
 Eigen::Array<double, N, 1> *u_errbd = 0)
{
    reorder_diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    u.transposeInPlace();
}

/**
 * Diagonalizes N-by-N symmetric matrix m so that
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * u
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 * @param[out] u      N-by-N complex unitary matrix
 */
template<class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u)
{
    fs_diagonalize_symmetric_errbd(m, s, u);
}

/**
 * Same as fs_diagonalize_symmetric(m, s, u) except that an
 * approximate error bound for the singular values is returned.  The
 * error bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd)
{
    fs_diagonalize_symmetric_errbd(m, s, u, &s_errbd);
}

/**
 * Same as fs_diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of fs_diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& s,
 Eigen::Matrix<std::complex<double>, N, N>& u,
 double& s_errbd,
 Eigen::Array<double, N, 1>& u_errbd)
{
    fs_diagonalize_symmetric_errbd(m, s, u, &s_errbd, &u_errbd);
}

template<int N>
struct CompareAbs {
    CompareAbs(const Eigen::Array<double, N, 1>& w_) : w(w_) {}
    bool operator() (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); }
    const Eigen::Array<double, N, 1>& w;
};

template<class Scalar, int N>
void fs_diagonalize_hermitian_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double *w_errbd = 0,
 Eigen::Array<double, N, 1> *z_errbd = 0)
{
    diagonalize_hermitian_errbd(m, w, z, w_errbd, z_errbd);
    Eigen::PermutationMatrix<N> p;
    p.setIdentity();
    std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
	      CompareAbs<N>(w));
#if EIGEN_VERSION_AT_LEAST(3,1,4)
    w.matrix().transpose() *= p;
    if (z_errbd) z_errbd->matrix().transpose() *= p;
#else
    Eigen::Map<Eigen::Matrix<double, N, 1> >(w.data()).transpose() *= p;
    if (z_errbd)
	Eigen::Map<Eigen::Matrix<double, N, 1> >(z_errbd->data()).transpose()
	    *= p;
#endif
    z = (z * p).adjoint().eval();
}

/**
 * Diagonalizes N-by-N hermitian matrix m so that
 *
 *     m == z.adjoint() * w.matrix().asDiagonal() * z    // convention of SARAH
 *
 * w is arranged so that `abs(w[i])` are in ascending order.
 *
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 * @param[out] z      N-by-N unitary matrix
 */
template<class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z)
{
    fs_diagonalize_hermitian_errbd(m, w, z);
}

/**
 * Same as fs_diagonalize_hermitian(m, w, z) except that an
 * approximate error bound for the eigenvalues is returned.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of fs_diagonalize_hermitian(m, w, z) for the
 * other parameters.
 */
template<class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double& w_errbd)
{
    fs_diagonalize_hermitian_errbd(m, w, z, &w_errbd);
}

/**
 * Same as fs_diagonalize_hermitian(m, w, z, w_errbd) except that
 * approximate error bounds for the eigenvectors are returned.  The
 * error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] z_errbd array of approximate error bounds for z
 *
 * See the documentation of fs_diagonalize_hermitian(m, w, z, w_errbd)
 * for the other parameters.
 */
template<class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<double, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 double& w_errbd,
 Eigen::Array<double, N, 1>& z_errbd)
{
    fs_diagonalize_hermitian_errbd(m, w, z, &w_errbd, &z_errbd);
}

} // namespace flexiblesusy

#endif // linalg2_hpp
