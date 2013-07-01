
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_betafunction

#include <boost/test/unit_test.hpp>

#include "rge.h"
#include "betafunction.hpp"

Eigen::ArrayXd ToEigenArray(const DoubleVector& v)
{
   Eigen::ArrayXd a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a(i - 1) = v(i);
   return a;
}

DoubleVector ToDoubleVector(const Eigen::ArrayXd& a)
{
   DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

class RGE_model : public RGE {
public:
   RGE_model() : RGE(), pars(10) {
      setPars(10);
      setMu(100.);
   }
   virtual ~RGE_model() {}
   virtual const DoubleVector display() const { return pars; }
   virtual void set(const DoubleVector& s) { pars = s; }
   virtual DoubleVector beta() const {
      DoubleVector beta(pars.size());
      for (std::size_t i = 1; i <= pars.size(); i++)
         beta(i) = pars(i) * pars(i) * 0.1 + 0.05 * pars(pars.size() + 1 - i);
      return beta;
   }
private:
   DoubleVector pars;
};

class Eigen_model : public Beta_function {
public:
   Eigen_model() : Beta_function(), pars(Eigen::ArrayXd::Zero(10)) {
      set_scale(100.);
      set_parameters(10);
   }
   virtual ~Eigen_model() {}
   virtual const Eigen::ArrayXd display() const { return pars; }
   virtual void set(const Eigen::ArrayXd& s) { pars = s; }
   virtual Eigen::ArrayXd beta() const {
      Eigen::ArrayXd beta(pars.size());
      for (int i = 0; i < pars.rows(); i++)
         beta(i) = pars(i) * pars(i) * 0.1 + 0.05 * pars(pars.rows() - 1 - i);
      return beta;
   }
private:
   Eigen::ArrayXd pars;
};

BOOST_AUTO_TEST_CASE( test_running )
{
   DoubleVector pars(10);
   for (std::size_t i = 1; i <= pars.size(); i++)
      pars(i) = 0.1 * i;

   RGE_model rge;
   Eigen_model eig;

   rge.set(pars);
   eig.set(ToEigenArray(pars));

   BOOST_CHECK_EQUAL(rge.display(), ToDoubleVector(eig.display()));
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
   BOOST_CHECK_EQUAL(rge.howMany(), eig.get_parameters());

   rge.runto(1.0e10);
   eig.run_to(1.0e10);

   BOOST_CHECK_EQUAL(rge.display(), ToDoubleVector(eig.display()));
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
}
