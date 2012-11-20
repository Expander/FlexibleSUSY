
#ifndef GUT_SCALE_RUNNER_H
#define GUT_SCALE_RUNNER_H

#include <cmath>

/**
 * @class GUT_scale_runner
 * @brief calculates the approximate GUT scale
 *
 * In the default implementation the GUT scale is defined by the scale
 * where g1 == g2.  The closer the current scale is to the GUT scale,
 * the better the guess is.
 *
 * Usage:
 *
 * @code
 * StandardModel* sm = new StandardModel();
 * GUT_scale_calculator<StandardModel> gsr;
 * double MX = 1.0e15; // guessed GUT scale
 * while (...) {
 *    sm->runto(MX);
 *    // recalculate GUT scale
 *    MX = gsr.calculateGUTScale(*sm);
 *    // ...
 * }
 * @endcode
 */
template <typename T>
class GUT_scale_calculator {
public:
   double calculateGUTScale(const T&) const;
};

template <typename T>
double GUT_scale_calculator<T>::calculateGUTScale(const T& rge) const
{
   const double currentScale = rge.getScale();
   const double g1 = rge.displayGaugeCoupling(1);
   const double g2 = rge.displayGaugeCoupling(2);
   const T beta(rge.calcBeta());
   const double betaG1 = beta.displayGaugeCoupling(1);
   const double betaG2 = beta.displayGaugeCoupling(2);

   const double gutScale = currentScale * exp((g1 - g2)/(betaG2 - betaG1));

   return gutScale;
}

#endif
