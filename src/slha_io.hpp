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

#ifndef SLHA_IO_H
#define SLHA_IO_H

#include <string>
#include <functional>

namespace softsusy {
   class QedQcd;
}

namespace flexiblesusy {

class SLHA_io {
public:
   typedef std::function<void(int, double)> Tuple_processor;

   SLHA_io();
   SLHA_io(const std::string&);
   virtual ~SLHA_io() {}

   const std::string& get_input_file() const { return input_filename; }
   void set_input_file(const std::string& f) { input_filename = f; }

   void fill(softsusy::QedQcd&);
   void read_block(const std::string&, Tuple_processor);

private:
   std::string input_filename;
};

} // namespace flexiblesusy

#endif
