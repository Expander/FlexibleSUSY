
#ifndef TWO_SCALE_CONSTRAINT_H
#define TWO_SCALE_CONSTRAINT_H

class DoubleVector;
class Two_scale;

template<>
class Constraint<Two_scale> {
public:
   virtual void apply(DoubleVector&) = 0;
};

typedef Constraint<Two_scale> Two_scale_constraint;

#endif
