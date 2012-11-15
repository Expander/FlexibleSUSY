
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_logger

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "logger.hpp"

BOOST_AUTO_TEST_CASE( test_all_messages )
{
   VERBOSE_MSG("verbose message")
   DEBUG_MSG("debug message")
   INFO("info message")
   WARNING("warning message")
   ERROR("error message")
   FATAL("fatal message")
}

BOOST_AUTO_TEST_CASE( test_streaming_operator )
{
   INFO("streamed info" << " message")
}
