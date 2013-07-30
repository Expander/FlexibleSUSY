// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef ERROR_H
#define ERROR_H

#include <string>

namespace flexiblesusy {

class Error {
public:
   virtual ~Error() {}
   virtual std::string what() const = 0;
};

class Two_scale_model;

class TachyonError : public Error {
public:
   TachyonError(const Two_scale_model*, const std::string&, int);
   virtual ~TachyonError() {}
   virtual std::string what() const;
   const Two_scale_model* get_model() const { return model; }
   const std::string& get_particle_name() const { return particle_name; }
   int get_particle_index() const { return particle_index; }
private:
   const Two_scale_model* model;
   std::string particle_name;
   int particle_index;
};

class NoEWSBError : public Error {
public:
   NoEWSBError(const Two_scale_model*, double);
   virtual ~NoEWSBError() {}
   virtual std::string what() const;
   const Two_scale_model* get_model() const { return model; }
   double get_requested_precision() const { return requested_precision; }
private:
   const Two_scale_model* model;
   double requested_precision;
};

}

#endif
