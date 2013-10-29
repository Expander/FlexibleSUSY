
#ifndef TEST_H
#define TEST_H

#include "linalg.h"
#include <Eigen/Core>

static const double max_dev = 1.0e-12;
static int gErrors = 0;

template <typename T>
bool is_zero(T a)
{
   return std::fabs(a) < std::numeric_limits<T>::epsilon();
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   return std::fabs(a - b) < prec;
}

template <typename T>
bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   if (std::fabs(a) < std::numeric_limits<T>::epsilon())
      return std::numeric_limits<T>::infinity();

   return std::fabs((a - b)/a) < prec;
}

bool is_equal(const DoubleMatrix& a, const DoubleMatrix& b, double max_dev)
{
   if (a.displayRows() != b.displayRows()) {
      std::cout << "matrices have different number of rows\n";
      return false;
   }
   if (a.displayCols() != b.displayCols()) {
      std::cout << "matrices have different number of columns\n";
      return false;
   }
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         if (!is_equal(a(i,l), b(i,l), max_dev))
            return false;
      }
   }
   return true;
}

bool is_equal(const Eigen::Matrix<double,3,3>& a, const Eigen::Matrix<double,3,3>& b, double max_dev)
{
   for (int i = 0; i < 3; ++i) {
      for (int l = 0; l < 3; ++l) {
         if (!is_equal(a(i,l), b(i,l), max_dev))
            return false;
      }
   }
   return true;
}

template <typename T>
void check_equality(T a, T b, const std::string& testMsg, T max_dev)
{
   if (!is_equal(a, b, max_dev)) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b
                << " (diff.: " << (a-b) << ", rel. diff.: "
                << 100. * (a-b)/a << "%)\n";
      gErrors++;
   }
}

void check_equality(Complex a, Complex b, const std::string& testMsg, double max_dev)
{
   std::ostringstream msg;
   msg << testMsg << " (real part)";
   check_equality(a.real(), b.real(), msg.str(), max_dev);
   msg.str("");
   msg << testMsg << " (imag part)";
   check_equality(a.imag(), b.imag(), msg.str(), max_dev);
}

void check_equality(int a, int b, const std::string& testMsg, double)
{
   if (a != b) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b << ")\n";
      gErrors++;
   }
}

void check_equality(const DoubleVector& a, const DoubleVector& b,
                    const std::string& testMsg, double max_dev)
{
   if (a.displayStart() != b.displayStart()) {
      std::cout << "test failed: " << testMsg << ": "
                << "vectors have different start value: "
                << "a.displayStart() == " << a.displayStart()
                << ", b.displayStart() == " << b.displayStart()
                << std::endl;
      gErrors++;
      return;
   }

   if (a.displayEnd() != b.displayEnd()) {
      std::cout << "test failed: " << testMsg << ": "
                << "vectors have different end value: "
                << "a.displayEnd() == " << a.displayEnd()
                << ", b.displayEnd() == " << b.displayEnd()
                << std::endl;
      gErrors++;
      return;
   }

   for (int i = a.displayStart(); i <= a.displayEnd(); ++i) {
      std::ostringstream element;
      element << testMsg << " [element " << i << "]";
      check_equality(a(i), b(i), element.str(), max_dev);
   }
}

void check_equality(const DoubleMatrix& a, const DoubleMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const DoubleMatrix& b,
                    const Eigen::Matrix<double,3,3>& a,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 0; i < 3; ++i) {
      for (int l = 0; l < 3; ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i+1,l+1), element.str(), max_dev);
      }
   }
}

void check_equality(const Eigen::Matrix<double,3,3>& a,
                    const Eigen::Matrix<double,3,3>& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 0; i < 3; ++i) {
      for (int l = 0; l < 3; ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const ComplexMatrix& a, const ComplexMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const DoubleMatrix& a, const ComplexMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "] ";
         check_equality(Complex(a(i,l), 0.0), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const ComplexMatrix& a, const DoubleMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   check_equality(b, a, testMsg, max_dev);
}

void check_condition(bool condition, const std::string& testMsg)
{
   if (!condition) {
      std::cout << "test failed: " << testMsg << "\n";
      gErrors++;
   }
}

template <typename T>
void check_greater_than(T a, T b, const std::string& testMsg)
{
   if (!(a > b)) {
      std::cout << "test failed: " << testMsg << ": " << a << " > "
                << b << "\n";
      gErrors++;
   }
}

template <typename T>
void check_relative_dev(T a, T b, const std::string& testMsg, T max_dev)
{
   if (!is_equal_rel(a, b, max_dev)) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b
                << " (diff.: " << (a-b) << ", rel. diff.: "
                << 100. * (a-b)/a << "%)\n";
      gErrors++;
   }
}

#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)
#define TEST_EQUALITY(a, b) check_equality(a, b,"line " S__LINE__ ": " #a " == " #b, max_dev)
#define TEST_CLOSE(a, b, dev) check_equality(a, b, "line " S__LINE__ ": " #a " == " #b, dev)
#define TEST_GREATER(a, b) check_greater_than(a, b,"line " S__LINE__ ": " #a " > " #b)
#define TEST_CLOSE_REL(a, b, dev) check_relative_dev(a, b, "line " S__LINE__ ": " #a " =r= " #b, dev)
#define TEST(condition) check_condition(condition, #condition);

#endif
