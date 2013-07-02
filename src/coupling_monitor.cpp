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

#include "coupling_monitor.hpp"

#include <fstream>
#include <iostream>

Coupling_monitor::Coupling_monitor()
   : couplings(TData())
{
}

Coupling_monitor::~Coupling_monitor()
{
}

/**
 * Get the couplings at the largest scale
 *
 * @return a pair with the scale and a Eigen::ArrayXd which contains the
 * couplings at this scale
 */
Coupling_monitor::TTouple Coupling_monitor::get_max_scale() const
{
   if (couplings.empty()) {
      ERROR("Data container is empty!");
      return TTouple(0.0, Eigen::ArrayXd(1));
   }

   // find gauge couplings at the greatest scale
   TData::const_iterator maxScale
      = max_element(couplings.begin(), couplings.end(), TDataComp());

   return *maxScale;
}

/**
 * Delete all internal couplings.
 */
void Coupling_monitor::reset()
{
   couplings.clear();
}

/**
 * write help line which describes the written data
 *
 * @param coupling_name coupling name
 * @param fout output stream
 * @param number_of_couplings number of couplings at each scale
 * @param start_index index of first coupling (usually set to 1)
 */
void Coupling_monitor::write_comment_line(char coupling_name, std::ofstream& fout,
                                          std::size_t number_of_couplings,
                                          int start_index) const
{
   if (!fout.good())
      return;

   fout << "# scale [GeV]";

   for (std::size_t i = 0; i < number_of_couplings; ++i)
      fout << " | " << coupling_name << "(" << start_index + i << ")";

   fout << std::endl;
}

/**
 * Write all couplings to a text file.
 *
 * @param file_name name of file to write the data to
 */
void Coupling_monitor::write_to_file(const std::string& file_name) const
{
   const char coupling_name = 'g';

   if (couplings.empty())
      return;

   std::ofstream filestr(file_name.c_str(), std::ios::out);
   VERBOSE_MSG("opening file: " << file_name.c_str());
   if (filestr.fail()) {
      ERROR("can't open file " << file_name
            << " for writing running couplings");
      return;
   }

   write_comment_line(coupling_name, filestr,
                      couplings.front().second.size(), 0);

   // write data
   for (TData::const_iterator it = couplings.begin();
        it != couplings.end(); ++it) {
      if (!filestr.good()) {
         ERROR("file " << file_name << " is corrupted");
         break;
      }

      filestr.width(16);
      filestr << std::left << it->first;

      // write all gauge couplings in order
      for (int i = 0; i < it->second.size(); ++i) {
         filestr.width(16);
         filestr << std::left << it->second(i);
      }

      filestr << std::endl;
   }

   filestr.close();
   VERBOSE_MSG("file written: " << file_name.c_str());
}
