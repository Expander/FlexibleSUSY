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
#include "eigen_utils.hpp"

using namespace flexiblesusy;

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

   Eigen::Array<double,3,1> dst;

   remove_if_equal(src, cmp, dst);

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

   Eigen::Array<double,3,1> dst;

   remove_if_equal(src, cmp, dst);

   BOOST_CHECK_EQUAL(dst(0), 1);
   BOOST_CHECK_EQUAL(dst(1), 4);
   BOOST_CHECK_EQUAL(dst(2), 5);
}
