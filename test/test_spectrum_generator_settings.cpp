#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_spectrum_generator_settings

#include <boost/test/unit_test.hpp>
#include "spectrum_generator_settings.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_get )
{
   Spectrum_generator_settings s;
   const auto es = s.get();

   BOOST_REQUIRE_EQUAL(Spectrum_generator_settings::NUMBER_OF_OPTIONS, es.size());

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      BOOST_CHECK_EQUAL(s.get(static_cast<Spectrum_generator_settings::Settings>(i)), es(i));
   }
}

BOOST_AUTO_TEST_CASE( test_set )
{
   const int N = Spectrum_generator_settings::NUMBER_OF_OPTIONS;
   Spectrum_generator_settings::Settings_t es = Spectrum_generator_settings::Settings_t::Random(N);

   BOOST_REQUIRE_EQUAL(Spectrum_generator_settings::NUMBER_OF_OPTIONS, es.size());

   Spectrum_generator_settings s;
   s.set(es);

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      BOOST_CHECK_EQUAL(s.get(static_cast<Spectrum_generator_settings::Settings>(i)), es(i));
   }
}
