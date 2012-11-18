
#ifndef TWO_SCALE_MATCHING_H
#define TWO_SCALE_MATCHING_H

#include "matching.hpp"

class Two_scale;

template<>
class Matching<Two_scale> {
public:
   virtual ~Matching() {}
   virtual void matchLowToHighScaleModel() const = 0;
   virtual void matchHighToLowScaleModel() const = 0;
};

#endif
