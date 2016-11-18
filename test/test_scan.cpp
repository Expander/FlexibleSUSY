#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_scan

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "scan.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( uniform_distribution_0_1 )
{
   const std::size_t N = 1000;
   std::array<double, N> a{};

   std::generate(std::begin(a), std::end(a),
                 [] () { return Uniform<>::dice(); });

   BOOST_CHECK(*std::max_element(std::begin(a), std::end(a)) <= 1.);
   BOOST_CHECK(*std::min_element(std::begin(a), std::end(a)) >= 0.);
}

BOOST_AUTO_TEST_CASE( uniform_distribution_start_stop )
{
   const std::size_t N = 1000;
   const double start = 0., stop = 10.;
   std::array<double, N> a{};

   std::generate(std::begin(a), std::end(a),
                 [start, stop] () { return Uniform<>::dice(start, stop); });

   BOOST_CHECK(*std::max_element(std::begin(a), std::end(a)) <= stop);
   BOOST_CHECK(*std::min_element(std::begin(a), std::end(a)) >= start);
}
