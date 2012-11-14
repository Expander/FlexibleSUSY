
#ifndef TWO_SCALE_MATCHING_H
#define TWO_SCALE_MATCHING_H

#include "matching.hpp"

class DoubleVector;
class Two_scale;

template<>
class Matching<Two_scale> {
public:
   virtual ~Matching() {}
   virtual void matchLowToHighScaleModel() const = 0;
   virtual void matchHighToLowScaleModel() const = 0;
};

typedef Matching<Two_scale> Two_scale_matching;

#endif
