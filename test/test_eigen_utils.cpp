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
#define BOOST_TEST_MODULE test_eigen_utils

#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include "eigen_utils.hpp"
#include <complex>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_div_safe_zero)
{
   Eigen::Array<double,3,3> v1, v2, result;
   v1.setZero();
   v2.setZero();
   result = div_safe(v1, v2);

   BOOST_CHECK(result.allFinite());

   for (int i = 0; i < result.rows(); i++)
      for (int k = 0; k < result.cols(); k++)
         BOOST_CHECK_EQUAL(result(i,k), 0.);
}

BOOST_AUTO_TEST_CASE(test_div_safe)
{
   Eigen::Array<double,3,1> v1, v2, result;
   v1 << 1, 2, 3;
   v2 << 0, 1, 2;
   result = div_safe(v1, v2);

   BOOST_CHECK_EQUAL(result(0), 0.);
   BOOST_CHECK_EQUAL(result(1), 2.);
   BOOST_CHECK_EQUAL(result(2), 3./2.);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector)
{
   Eigen::Array<double,3,1> vec1, vec2;
   vec1 << 1, 2, 3;
   vec2 << 3, 1, 2;

   reorder_vector(vec1, vec2);

   BOOST_CHECK_EQUAL(vec1(0), 3);
   BOOST_CHECK_EQUAL(vec1(1), 1);
   BOOST_CHECK_EQUAL(vec1(2), 2);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector_abs)
{
   Eigen::Array<double,3,1> vec1, vec2;
   vec1 << 1, 2, 3;
   vec2 << 1, 2, -3;

   reorder_vector(vec1, vec2);

   BOOST_CHECK_EQUAL(vec1(0), 2);
   BOOST_CHECK_EQUAL(vec1(1), 3);
   BOOST_CHECK_EQUAL(vec1(2), 1);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector_with_matrix)
{
   Eigen::Array<double,3,1> vec1;
   vec1 << 1, 2, 3;

   Eigen::Matrix<double,3,3> matrix;
   matrix << 3, 0, 0,
             0, 1, 0,
             0, 0, 2;

   reorder_vector(vec1, matrix);

   BOOST_CHECK_EQUAL(vec1(0), 3);
   BOOST_CHECK_EQUAL(vec1(1), 1);
   BOOST_CHECK_EQUAL(vec1(2), 2);
}

BOOST_AUTO_TEST_CASE(test_remove_if_equal)
{
   Eigen::Array<double,5,1> src;
   src << 1, 2, 3, 4, 5;

   Eigen::Array<double,2,1> cmp; // to be removed from src
   cmp << 2, 4;

   const Eigen::Array<double,3,1> dst = remove_if_equal(src, cmp);

   BOOST_CHECK_EQUAL(dst(0), 1);
   BOOST_CHECK_EQUAL(dst(1), 3);
   BOOST_CHECK_EQUAL(dst(2), 5);
}

BOOST_AUTO_TEST_CASE(test_remove_if_equal_double_indices)
{
   Eigen::Array<double,5,1> src;
   src << 1, 2, 2.3, 4, 5;

   Eigen::Array<double,2,1> cmp; // to be removed from src
   cmp << 2, 2.1;

   const Eigen::Array<double,3,1> dst = remove_if_equal(src, cmp);

   BOOST_CHECK_EQUAL(dst(0), 1);
   BOOST_CHECK_EQUAL(dst(1), 4);
   BOOST_CHECK_EQUAL(dst(2), 5);
}

BOOST_AUTO_TEST_CASE(test_Eval)
{
   Eigen::Matrix<double,2,2> m1, m2;
   Eigen::Matrix<double,2,1> v1, v2;
   int i;
   double d;
   std::complex<double> c;

   BOOST_CHECK_EQUAL(typeid(Eval(m1)).hash_code(), typeid(m1).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(m1*m2)).hash_code(), typeid(m1).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(m1*m2 + m1)).hash_code(), typeid(m1).hash_code());

   BOOST_CHECK_EQUAL(typeid(Eval(v1)).hash_code(), typeid(v1).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(m1*v1)).hash_code(), typeid(v1).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(m1*v1 + v2)).hash_code(), typeid(v1).hash_code());

   BOOST_CHECK_EQUAL(typeid(Eval(c)).hash_code(), typeid(c).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(d)).hash_code(), typeid(d).hash_code());
   BOOST_CHECK_EQUAL(typeid(Eval(i)).hash_code(), typeid(i).hash_code());

}
