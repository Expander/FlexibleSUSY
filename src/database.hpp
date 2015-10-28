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

#include <string>
#include <Eigen/Core>

struct sqlite3;

namespace flexiblesusy {

class Database {
public:
   Database(const std::string& file_name);
   ~Database();

   void insert(const std::string&, const std::vector<std::string>&, const Eigen::ArrayXd&);

private:
   sqlite3* db;

   sqlite3* open(const std::string&);
   void execute(const std::string&);
   template <typename T> void create_table(const std::string&, const std::vector<std::string>&);
};

} // namespace flexiblesusy
