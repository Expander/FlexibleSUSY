
#include "wrappers.hpp"

double AbsSqr(const Complex& z)
{
   return z.norm();
}

double ArcTan(double a)
{
   return atan(a);
}

double ArcSin(double a)
{
   return asin(a);
}

double ArcCos(double a)
{
   return acos(a);
}

DoubleMatrix Adj(const DoubleMatrix& m)
{
   return m.transpose();
}

ComplexMatrix Adj(const ComplexMatrix& m)
{
   return m.hermitianConjugate();
}

double Conj(double a)
{
   return a;
}

Complex Conj(const Complex& a)
{
   return a.conj();
}

DoubleMatrix Conj(const DoubleMatrix& m)
{
   return m;
}

ComplexMatrix Conj(const ComplexMatrix& m)
{
   return m.complexConjugate();
}

double Cos(double x)
{
   return cos(x);
}

double Sin(double x)
{
   return sin(x);
}

int Delta(int i, int j)
{
   return i == j ? 1 : 0;
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
 * diag = u m u^+
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

/**
 * diag = u m u^T
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (unsorted)
 */
void DiagonalizeUnsorted(const DoubleMatrix& m, DoubleMatrix& u,
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
   u = u.transpose();
}

/**
 * diag = u m u^+
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (unsorted)
 */
void DiagonalizeUnsorted(const DoubleMatrix& m, ComplexMatrix& u,
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
   DoubleMatrix a(m), k(c, c);
   int nrot;
   diagonaliseJac(a, c, eigenvalues, k, &nrot);

   // We want to change the PHASES of the neutralino mixing matrix in order to
   // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
   const ComplexMatrix trans(phase_rotation(eigenvalues));
   u = trans * k.transpose();
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const DoubleMatrix& m, DoubleMatrix& u , DoubleVector& eigenvalues)
{
   double theta;
   eigenvalues = m.sym2by2(theta);
   u = rot2d(theta);
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const DoubleMatrix& m, ComplexMatrix& u, DoubleVector& eigenvalues)
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

/**
 * diag = u m v^T
 *
 * @param m mass matrix
 * @param u mixing matrix
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

   // Phase freedom - tends to give more familiar matrices
   u = -1.0 * u; v = -1.0 * v;

   AssociateOrderAbs(u, v, eigenvalues);

   v = v.transpose();
   u = u.transpose();
}

/**
 * diag = u m v^+
 *
 * @param m mass matrix
 * @param u mixing matrix
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

double Log(double a)
{
   return std::log(a);
}

double Mass2(double m)
{
   return m * m;
}

double MaxRelDiff(double a, double b)
{
   return toleranceCheck(a, b);
}

double MaxRelDiff(const DoubleVector& a, const DoubleVector& b)
{
   return a.compare(b);
}

double Re(double x)
{
   return x;
}

double Re(const Complex& x)
{
   return real(x);
}

int ThetaStep(int a, int b)
{
   return a <= b ? 1 : 0;
}

DoubleMatrix Tp(const DoubleMatrix& m)
{
   return m.transpose();
}

ComplexMatrix Tp(const ComplexMatrix& m)
{
   return m.transpose();
}

DoubleMatrix Transpose(const DoubleMatrix& m)
{
   return m.transpose();
}

ComplexMatrix Transpose(const ComplexMatrix& m)
{
   return m.transpose();
}

double trace(const DoubleMatrix& m)
{
   return m.trace();
}

Complex trace(const ComplexMatrix& m)
{
   return m.trace();
}

double Sqrt(double a)
{
   return std::sqrt(a);
}
