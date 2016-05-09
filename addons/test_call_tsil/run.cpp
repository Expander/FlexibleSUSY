#include "call_tsil.h"
#include "call_tsil.hpp"
#include <iostream>

int main()
{
   std::cout << "C  : A0(100,100) = " << call_A(100., 100.) << '\n';
   std::cout << "C++: A0(100,100) = " << call_A_cpp(100., 100.) << '\n';

   return 0;
}
