
#ifndef TWO_SCALE_MATCHING_H
#define TWO_SCALE_MATCHING_H

#include "matching.hpp"

namespace flexiblesusy {

class Model;
class Two_scale;

template<>
class Matching<Two_scale> {
public:
   virtual ~Matching() {}
   virtual void match() = 0;
   virtual double get_scale() const = 0;
   virtual void set_models(Model*, Model*) = 0;
};

}

#endif
