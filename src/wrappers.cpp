
#include "wrappers.hpp"

#include <Eigen/SVD>

namespace flexiblesusy {

Eigen::Matrix3d Diag(const Eigen::Matrix3d& m)
{
   Eigen::Matrix3d diag(m);
   for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
         if (i != k)
            diag(i,k) = 0.0;
      }
   }
   return diag;
}

/**
 * diag = u m u^T
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const DoubleMatrix& m, DoubleMatrix& u,
                 DoubleVector& eigenvalues)
{
   const int c = m.displayCols();
#ifdef DEBUG
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       u.displayCols() != c || u.displayRows() !=c || c > 10) {
      ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw ii.str();
   }
#endif
   // Numerical recipes routine replaces argument, so make a copy of elements
   DoubleMatrix a(m);
   int nrot;
   diagonaliseJac(a, c, eigenvalues, u, &nrot);
   u.associateOrderAbs(eigenvalues);
   u = u.transpose();
}

/**
 * Returns a phase rotation matrix which will transform a mixing
 * matrix such that the corresponding eigenvalues are positive.
 *
 * @param v eigen values
 *
 * @return phase rotation matrix
 */
ComplexMatrix phase_rotation(const DoubleVector& v)
{
   const int len = v.size();
   ComplexMatrix rot(len, len);

   for (int i = 1; i <= len; i++) {
      if (v(i) < 0.0)
         rot(i, i) = Complex(0.0, -1.0);
      else
         rot(i, i) = Complex(1.0, 0.0);
   }

   return rot;
}

/**
 * diag = u^* m u^+
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const DoubleMatrix& m, ComplexMatrix& u, DoubleVector& eigenvalues)
{
   const int c = m.displayCols();
#ifdef DEBUG
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       u.displayCols() != c || u.displayRows() !=c || c > 10) {
      ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw ii.str();
   }
#endif
   // Numerical recipes routine replaces argument, so make a copy of elements
   DoubleMatrix a(m), k(c, c);
   int nrot;
   diagonaliseJac(a, c, eigenvalues, k, &nrot);

   // So Far, V is such that the eigenvalues are in some random order.
   // We must now re-order the rows of V to get the eigenvalues in the
   // correct order (taken to be increasing abs(eigenvalues)).
   k.associateOrderAbs(eigenvalues);

   // We want to change the PHASES of the neutralino mixing matrix in order to
   // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
   const ComplexMatrix trans(phase_rotation(eigenvalues));
   u = trans * k.transpose();
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const DoubleMatrix& m, DoubleMatrix& u,
                     DoubleVector& eigenvalues)
{
   double theta;
   eigenvalues = m.sym2by2(theta);
   if (std::fabs(eigenvalues(1)) > std::fabs(eigenvalues(2))) {
      theta += 0.5 * M_PI;
      std::swap(eigenvalues(1), eigenvalues(2));
   }
   u = rot2d(theta);
}

void Diagonalize2by2(const DoubleMatrix& m, ComplexMatrix& u,
                     DoubleVector& eigenvalues)
{
   double theta;
   eigenvalues = m.sym2by2(theta);
   const ComplexMatrix a(phase_rotation(eigenvalues).complexConjugate());
   u = a * rot2d(theta);
   eigenvalues = eigenvalues.apply(fabs);
}

void AssociateOrderAbs(DoubleMatrix& u, DoubleMatrix& v, DoubleVector& w)
{
#ifdef DEBUG
   if ((u.displayRows() != w.displayEnd()) ||
       (u.displayCols() != w.displayEnd()) ||
       (w.displayStart() != 1)) {
      ostringstream ii;
      ii << "Associated ordering incompatibility\n";
      throw ii.str();
   }
#endif
   for (int i = w.displayStart(); i <= w.displayEnd(); ++i) {
      for (int j = i + 1; j <= w.displayEnd(); ++j) {
         if (abs(w.display(i)) > abs(w.display(j))) {
            w.swap(i, j);
            v.swapcols(i, j);
            u.swapcols(i, j);
         }
      }
   }
}

Eigen::MatrixXd ToEigen(const DoubleMatrix& m)
{
   const int cols = m.displayCols(), rows = m.displayRows();
   Eigen::MatrixXd eig(rows,cols);
   for (int i = 1; i <= rows; ++i)
      for (int j = 1; j <= rows; ++j)
         eig(i-1,j-1) = m(i,j);
   return eig;
}

DoubleMatrix ToSoftsusy(const Eigen::MatrixXd& m)
{
   const int cols = m.cols(), rows = m.rows();
   DoubleMatrix soft(rows,cols);
   for (int i = 1; i <= rows; ++i)
      for (int j = 1; j <= rows; ++j)
         soft(i,j) = m(i-1,j-1);
   return soft;
}

DoubleVector ToSoftsusy(const Eigen::VectorXd& m)
{
   const int cols = m.rows();
   DoubleVector soft(cols);
   for (int j = 1; j <= cols; ++j)
      soft(j) = m(j-1);
   return soft;
}

/**
 * diag = u m v^T
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param v mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const DoubleMatrix& m, DoubleMatrix& u,
                 DoubleMatrix& v, DoubleVector& eigenvalues)
{
#ifdef DEBUG
   const int c = m.displayCols();
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       v.displayCols() != c || v.displayRows() !=c ||
       u.displayCols() !=c || u.displayRows() !=c) {
      ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m
         << "with u" << u << "v " << v << "eigenvalues " << eigenvalues;
      throw ii.str();
   }
#endif
   // Numerical routine replaces argument, so make a copy of elements
   u = m;
   diagonaliseSvd(u, eigenvalues, v);

   // Eigen::MatrixXd me(ToEigen(m)), ue(ToEigen(u)), ve(ToEigen(v));
   // Eigen::JacobiSVD<Eigen::MatrixXd> eigensolver(me, Eigen::ComputeThinU | Eigen::ComputeThinV);
   // u = ToSoftsusy(eigensolver.matrixU());
   // v = ToSoftsusy(eigensolver.matrixV());
   // eigenvalues = ToSoftsusy(eigensolver.singularValues());

   // Phase freedom - tends to give more familiar matrices
   u = -1.0 * u; v = -1.0 * v;

   AssociateOrderAbs(u, v, eigenvalues);

   v = v.transpose();
   u = u.transpose();
}

/**
 * diag = u^* m v^+
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param v mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const DoubleMatrix& m, ComplexMatrix& u,
                 ComplexMatrix& v, DoubleVector& eigenvalues)
{
#ifdef DEBUG
   const int c = m.displayCols();
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       v.displayCols() != c || v.displayRows() !=c || c > 10) {
      ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw ii.str();
   }
#endif
   DoubleMatrix realU(u.real()), realV(v.real());
   Diagonalize(m, realU, realV, eigenvalues);

   const ComplexMatrix trans(phase_rotation(eigenvalues));
   u = trans * realU;
   v = trans * realV;
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const DoubleMatrix& m, ComplexMatrix& u,
                     ComplexMatrix& v, DoubleVector& eigenvalues)
{
   double thetaL, thetaR;
   eigenvalues = m.asy2by2(thetaL, thetaR);
   positivise(thetaL, thetaR, eigenvalues, u, v);
   eigenvalues = eigenvalues.apply(fabs);
}

double MaxRelDiff(double a, double b)
{
   const double sTin = fabs(a), sTout = fabs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return fabs(1.0 - std::min(sTin, sTout) / maxx);
}

double MaxRelDiff(const DoubleVector& a, const DoubleVector& b)
{
   return a.compare(b);
}

void Symmetrize(DoubleMatrix& m)
{
   const int r = m.displayRows();
   const int c = m.displayCols();
   for (int i = 1; i <= r; i++)
      for (int k = 1; k < i && k <= c; k++)
         m(i,k) = m(k,i);
}

Eigen::ArrayXd ToEigenArray(const DoubleVector& v)
{
   Eigen::ArrayXd a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a(i - v.displayStart()) = v(i);
   return a;
}

Eigen::ArrayXd ToEigenArray(double v)
{
   Eigen::ArrayXd a(1);
   a(0) = v;
   return a;
}

std::valarray<double> ToValarray(const DoubleVector& v)
{
   std::valarray<double> a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a[i - v.displayStart()] = v(i);
   return a;
}

std::valarray<double> ToValarray(double v)
{
   std::valarray<double> a(1);
   a[0] = v;
   return a;
}

DoubleVector ToDoubleVector(const Eigen::ArrayXd& a)
{
   DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

Eigen::MatrixXd ToEigenMatrix(const DoubleMatrix& m)
{
   const int r = m.displayRows();
   const int c = m.displayCols();
   Eigen::MatrixXd result(r,c);

   for (int i = 1; i <= r; i++)
      for (int k = 1; k <= c; k++)
         result(i-1, k-1) = m(i,k);

   return result;
}

DoubleMatrix ToDoubleMatrix(const Eigen::MatrixXd& m)
{
   const int r = m.rows();
   const int c = m.cols();
   DoubleMatrix result(r,c);

   for (int i = 0; i < r; i++)
      for (int k = 0; k < c; k++)
         result(i+1, k+1) = m(i,k);

   return result;
}

}
