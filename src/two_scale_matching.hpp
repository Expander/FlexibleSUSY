
#ifndef TWO_SCALE_MATCHING_H
#define TWO_SCALE_MATCHING_H

#include "matching.hpp"

namespace flexiblesusy {

class Two_scale;
class Two_scale_model;

template<>
class Matching<Two_scale> {
public:
   virtual ~Matching() = default;
   virtual void match() = 0;
   virtual double get_scale() const = 0;
   virtual void set_models(Two_scale_model*, Two_scale_model*) = 0;
};

}

#endif
