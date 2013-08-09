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
#include "slhaea.h"

namespace softsusy {
   class QedQcd;
}

namespace flexiblesusy {

class SLHA_io {
public:
   typedef std::function<void(int, double)> Tuple_processor;

   SLHA_io();
   SLHA_io(const std::string&);
   ~SLHA_io() {}

   const std::string& get_input_file() const { return input_filename; }
   void set_input_file(const std::string&);

   void fill(softsusy::QedQcd&);
   void read_block(const std::string&, Tuple_processor);
   void write_to_file(const std::string&);

private:
   std::string input_filename; ///< SHLA input file name
   SLHAea::Coll data;          ///< SHLA data
   static void process_sminputs_tuple(softsusy::QedQcd&, int, double);
};

} // namespace flexiblesusy

#endif
