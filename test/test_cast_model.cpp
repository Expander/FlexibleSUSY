
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cast_model

#include <boost/test/unit_test.hpp>

#include "two_scale_model.hpp"
#include "two_scale_constraint.hpp"

using namespace flexiblesusy;

class MyModel : public Two_scale_model {
public:
   virtual ~MyModel() {}
   virtual void calculate_spectrum() {}
   virtual std::string name() const { return "MyModel"; }
   virtual void run_to(double, double eps = -1.0) {}
   virtual void set_precision(double) {}
};

BOOST_AUTO_TEST_CASE( test_cast_model )
{
   MyModel model;

   Two_scale_model* model_ptr = &model;

   MyModel* mymodel_ptr = cast_model<MyModel*>(model_ptr);

   BOOST_CHECK(mymodel_ptr != NULL);
}
