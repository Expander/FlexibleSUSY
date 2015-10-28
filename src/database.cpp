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

#include "database.hpp"
#include "logger.hpp"

#include <limits>
#include <sstream>
#include <iomanip>
#include <sqlite3.h>

namespace flexiblesusy {

Database::Database(const std::string& file_name)
   : db(open(file_name))
{
}

Database::~Database()
{
   if (db)
      sqlite3_close(db);
}

void Database::insert(
   const std::string& table_name, const std::vector<std::string>& names,
   const Eigen::ArrayXd& data)
{
   const std::size_t number_of_elements = data.rows();

   if (names.size() != number_of_elements) {
      ERROR("number of names does not match vector size!");
      return;
   }

   create_table<double>(table_name, names);

   std::string sql("INSERT INTO " + table_name + " (");

   for (std::size_t i = 0; i < number_of_elements; i++) {
      sql += '"' + names[i] + '"';
      if (i + 1 != number_of_elements)
         sql += ',';
   }

   sql += ") VALUES (";

   for (std::size_t i = 0; i < number_of_elements; i++) {
      sql += to_string(data[i]);
      if (i + 1 != number_of_elements)
         sql += ',';
   }

   sql += ");";

   execute(sql);
}

template <typename T>
void Database::create_table(
   const std::string& table_name, const std::vector<std::string>& names)
{
   const std::size_t number_of_elements = names.size();
   std::string sql("CREATE TABLE IF NOT EXISTS " + table_name + " (");

   for (std::size_t i = 0; i < number_of_elements; i++) {
      sql += '"' + names[i] + '"' + " REAL";
      if (i + 1 != number_of_elements)
         sql += ',';
   }

   sql += ");";

   execute(sql);
}

void Database::execute(const std::string& cmd)
{
   char* zErrMsg = 0;
   const int rc = sqlite3_exec(db, cmd.c_str(), 0, 0, &zErrMsg);

   if (rc != SQLITE_OK) {
      ERROR("SQL error while executing command \"" << cmd << "\": " << zErrMsg);
      sqlite3_free(zErrMsg);
   } else {
      VERBOSE_MSG("SQL command \"" << cmd << "\" executed successfully");
   }
}

sqlite3* Database::open(const std::string& file_name)
{
   sqlite3* db = 0;
   char* zErrMsg = 0;

   const int rc = sqlite3_open(file_name.c_str(), &db);

   if (rc) {
      ERROR("Can't open database: " << sqlite3_errmsg(db));
      db = 0;
   }

   return db;
}

template <typename T>
std::string Database::to_string(T number)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<T>::digits10) << number;
    return out.str();
}

} // namespace flexiblesusy
