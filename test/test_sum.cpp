#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sum

#include <boost/test/unit_test.hpp>
#include "sum.hpp"
#include "stopwatch.hpp"
#include <complex>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_SUM_unsigned)
{
   BOOST_CHECK_EQUAL(SUM(i,1,0,1), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1), 1);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1), 2);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1), 3);

   BOOST_CHECK_EQUAL(SUM(i,1,0,i), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,i), 1);
   BOOST_CHECK_EQUAL(SUM(i,1,2,i), 3);
   BOOST_CHECK_EQUAL(SUM(i,1,3,i), 6);

   BOOST_CHECK_EQUAL(SUM(i,1,0,2*i+1), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,2*i+1), 3);
   BOOST_CHECK_EQUAL(SUM(i,1,2,2*i+1), 8);
   BOOST_CHECK_EQUAL(SUM(i,1,3,2*i+1), 15);

   BOOST_CHECK_EQUAL(SUM(i,1,2,SUM(k,1,2,k*i)), 9);
   BOOST_CHECK_EQUAL(SUM(i,1,2,i*SUM(k,1,2,k*i)), 15);
}

BOOST_AUTO_TEST_CASE(test_SUM_double)
{
   BOOST_CHECK_EQUAL(SUM(i,1,0,1.), 0.);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1.), 1.);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1.), 2.);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1.), 3.);

   BOOST_CHECK_EQUAL(SUM(i,1,0,1.*i), 0.);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1.*i), 1.);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1.*i), 3.);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1.*i), 6.);

   BOOST_CHECK_EQUAL(SUM(i,1,2,SUM(k,1,2,1.5*k*i)), 13.5);
}

BOOST_AUTO_TEST_CASE(test_SUM_complex)
{
   const std::complex<double> sum_1 = SUM(i,1,3,std::complex<double>(i,i));
   const std::complex<double> sum_2 = SUM(i,1,2,SUM(k,1,2,std::complex<double>(k*i,k*i)));

   BOOST_CHECK_EQUAL(sum_1.real(), 6.);
   BOOST_CHECK_EQUAL(sum_1.imag(), 6.);
   BOOST_CHECK_EQUAL(sum_2.real(), 9.);
   BOOST_CHECK_EQUAL(sum_2.imag(), 9.);
}

BOOST_AUTO_TEST_CASE(test_SUM_benchmark_1)
{
   const std::size_t N = 10000;
   Stopwatch stopwatch;

   stopwatch.start();
   const std::size_t result_SUM = SUM(i,1,N,SUM(k,1,N,k*i));
   stopwatch.stop();
   const double time_SUM = stopwatch.get_time_in_seconds();

   stopwatch.start();
   std::size_t result_loop = 0;
   for (std::size_t i = 1; i <= N; i++)
      for (std::size_t k = 1; k <= N; k++)
         result_loop += k*i;
   stopwatch.stop();
   const double time_loop = stopwatch.get_time_in_seconds();

   BOOST_MESSAGE("SUM  took " << time_SUM  << " seconds");
   BOOST_MESSAGE("loop took " << time_loop << " seconds");

   BOOST_CHECK_EQUAL(result_loop, result_SUM);

   // The passing of the following test justifies the expansion of the
   // sums of the Feynman diagrams at the meta code level, instead of
   // using the SUM() macro.  If this test fails at some later point,
   // we should use the SUM() macro, because it is much more readable.
   BOOST_CHECK_LT(10*time_loop, time_SUM);
}
