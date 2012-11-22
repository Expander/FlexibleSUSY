
#ifndef TWO_SCALE_MATCHING_H
#define TWO_SCALE_MATCHING_H

#include "matching.hpp"

class Two_scale;

template<>
class Matching<Two_scale> {
public:
   virtual ~Matching() {}
   virtual void match_low_to_high_scale_model() const = 0;
   virtual void match_high_to_low_scale_model() const = 0;
   virtual double get_scale() const = 0;
   virtual void update_scale() = 0;
};

#endif
