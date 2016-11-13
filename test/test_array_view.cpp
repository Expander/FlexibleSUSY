#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_array_view

#include <boost/test/unit_test.hpp>

#include "array_view.hpp"
#include <type_traits>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_empty_range )
{
   Dynamic_array_view<int> av;

   BOOST_CHECK_EQUAL(av.begin(), av.end());
   BOOST_CHECK_EQUAL(av.cbegin(), av.cend());
   BOOST_CHECK(av.data() == nullptr);
   BOOST_CHECK(av.empty());
   BOOST_CHECK_EQUAL(av.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_length_init )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, 4);

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_range_init )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, a + 4);

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_array_init )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a);

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_make_lenght_init )
{
   double a[4] = { 1., 2., 3., 4. };
   auto av = make_dynamic_array_view(a, 4);

   static_assert(std::is_same<decltype(av)::Element_t, double>(),
                 "Element type is not double");

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_make_range_init )
{
   double a[4] = { 1., 2., 3., 4. };
   auto av = make_dynamic_array_view(a, a + 4);

   static_assert(std::is_same<decltype(av)::Element_t, double>(),
                 "Element type is not double");

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_make_array_init )
{
   double a[4] = { 1., 2., 3., 4. };
   auto av = make_dynamic_array_view(a);

   static_assert(std::is_same<decltype(av)::Element_t, double>(),
                 "Element type is not double");

   BOOST_CHECK(!av.empty());
   BOOST_CHECK_EQUAL(av.size(), 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_readable )
{
   double a[4] = { 1., 2., 3., 4. };
   const Dynamic_array_view<double> av(a, a + 4);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_writeable )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, a + 4);

   av[0] += 1;
   av[1] += 1;
   av[2] += 1;
   av[3] += 1;

   BOOST_CHECK_EQUAL(av[0], 1. + 1.);
   BOOST_CHECK_EQUAL(av[1], 2. + 1.);
   BOOST_CHECK_EQUAL(av[2], 3. + 1.);
   BOOST_CHECK_EQUAL(av[3], 4. + 1.);

   BOOST_CHECK_EQUAL(av[0], a[0]);
   BOOST_CHECK_EQUAL(av[1], a[1]);
   BOOST_CHECK_EQUAL(av[2], a[2]);
   BOOST_CHECK_EQUAL(av[3], a[3]);
}

BOOST_AUTO_TEST_CASE( test_throw )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, a + 4);

   BOOST_CHECK_NO_THROW(av[0]);
   BOOST_CHECK_THROW(av[4], OutOfBoundsError);
   BOOST_CHECK_THROW(av[-1], OutOfBoundsError);
}

BOOST_AUTO_TEST_CASE( test_range_based_for )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, a + 4);

   for (auto x: av)
      BOOST_TEST_MESSAGE(x);
}

BOOST_AUTO_TEST_CASE( test_operator_output )
{
   double a[4] = { 1., 2., 3., 4. };
   Dynamic_array_view<double> av(a, a + 4);

   BOOST_TEST_MESSAGE(av);
}
