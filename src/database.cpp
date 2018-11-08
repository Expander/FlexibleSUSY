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
#include "config.h"

#include <limits>
#include <sstream>
#include <iomanip>

#ifdef ENABLE_SQLITE

#include <boost/lexical_cast.hpp>
#include <sqlite3.h>

namespace flexiblesusy {
namespace database {

Database::Database(const std::string& file_name)
   : db(open(file_name))
{
}

Database::~Database()
{
   sqlite3_close(db);
}

/**
 * Insert a row of doubles into a table.  If the table does not exist,
 * it is created.
 *
 * @param table_name name of table
 * @param names vector of column names
 * @param data vector of doubles
 */
void Database::insert(
   const std::string& table_name, const std::vector<std::string>& names,
   const Eigen::ArrayXd& data)
{
   const std::size_t number_of_elements = data.rows();

   if (names.size() != number_of_elements) {
      ERROR("number of column names (" << names.size() <<
            ") does not match vector size (" << number_of_elements << ")!");
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
      sql += boost::lexical_cast<std::string>(data[i]);
      if (i + 1 != number_of_elements)
         sql += ',';
   }

   sql += ");";

   execute("PRAGMA synchronous = OFF;");
   execute(sql);
}

/**
 * Extract a row of doubles from a table.
 *
 * @param table_name name of table
 * @param row row index (0 = 1st row, 1 = 2nd row, ..., -1 = last row,
 * -2 is 2nd to last row, ...)
 */
Eigen::ArrayXd Database::extract(const std::string& table_name, long long row)
{
   Eigen::ArrayXd values;
   const std::string sql =
      (row >= 0 ?
       "SELECT * FROM " + table_name + " LIMIT 1 OFFSET "
          + std::to_string(row) + ";" :
       "SELECT * FROM " + table_name + " WHERE ROWID = (SELECT MAX(ROWID) - "
          + std::to_string(std::abs(row + 1)) + " FROM " + table_name + ");");

   execute(sql, extract_callback, static_cast<void*>(&values));

   return values;
}

/**
 * Create a table in the database.  If the table already exists, this
 * function does nothing.
 *
 * @param table_name name of table
 * @param names names of table columns
 */
template <typename T>
void Database::create_table(
   const std::string& table_name, const std::vector<std::string>& names)
{
   const std::size_t number_of_elements = names.size();
   std::string sql("CREATE TABLE IF NOT EXISTS " + table_name + " (");

   for (std::size_t i = 0; i < number_of_elements; i++) {
      sql += '"' + names[i] + R"(" REAL)";
      if (i + 1 != number_of_elements)
         sql += ',';
   }

   sql += ");";

   execute(sql);
}

/**
 * Execute a SQL command without a callback.
 *
 * @param cmd command
 */
void Database::execute(const std::string& cmd)
{
   char* zErrMsg = nullptr;
   const int rc = sqlite3_exec(db, cmd.c_str(), nullptr, nullptr, &zErrMsg);

   if (rc != SQLITE_OK) {
      ERROR("SQL error while executing command \"" << cmd << "\": " << zErrMsg);
      sqlite3_free(zErrMsg);
   } else {
      VERBOSE_MSG("SQL command \"" << cmd << "\" executed successfully");
   }
}

/**
 * Execute a SQL command using a callback function.
 *
 * @param cmd command
 * @param callback pointer to callback function
 * @param data pointer to data passed to callback function
 */
void Database::execute(const std::string& cmd, TCallback callback, void* data)
{
   char* zErrMsg = nullptr;
   const int rc = sqlite3_exec(db, cmd.c_str(), callback, data, &zErrMsg);

   if (rc != SQLITE_OK) {
      ERROR("SQL error while executing command \"" << cmd << "\": " << zErrMsg);
      sqlite3_free(zErrMsg);
   } else {
      VERBOSE_MSG("SQL command \"" << cmd << "\" executed successfully");
   }
}

/**
 * Open database
 *
 * @param file_name file name
 *
 * @return pointer to database, or 0 if database cannot be opened.
 */
sqlite3* Database::open(const std::string& file_name)
{
   sqlite3* db = nullptr;

   const int rc = sqlite3_open(file_name.c_str(), &db);

   if (rc) {
      throw SQLiteReadError("Cannot open sqlite3 database file "
                            + file_name + ": " + sqlite3_errmsg(db));
   }

   return db;
}

/**
 * Callback function which fills an Eigen::ArrayXd with the data in
 * the given row.
 *
 * @param data pointer to Eigen::ArrayXd
 * @param argc number of columns
 * @param argv array of column entries
 * @param col_name array of column names
 *
 * @return 0
 */
#ifdef ENABLE_VERBOSE
int Database::extract_callback(void* data, int argc, char** argv, char** col_name)
#else
int Database::extract_callback(void* data, int argc, char** argv, char**)
#endif
{
   auto values = static_cast<Eigen::ArrayXd*>(data);
   values->conservativeResize(argc);

   for (int i = 0; i < argc; i++) {
      (*values)(i) = boost::lexical_cast<double>(argv[i]);
      VERBOSE_MSG(col_name[i] << " = " << argv[i]);
   }

   return 0;
}

} // namespace database
} // namespace flexiblesusy

#else

namespace flexiblesusy {
namespace database {

Database::Database(const std::string&)
   : db(nullptr)
{
}

Database::~Database()
{
}

void Database::insert(
   const std::string&, const std::vector<std::string>&, const Eigen::ArrayXd&)
{
   throw DisabledSQLiteError("Cannot call insert(), because SQLite support is disabled.");
}

Eigen::ArrayXd Database::extract(const std::string&, long long)
{
   throw DisabledSQLiteError("Cannot call extract(), because SQLite support is disabled.");
}

template <typename T>
void Database::create_table(const std::string&, const std::vector<std::string>&)
{
}

void Database::execute(const std::string&) {}

void Database::execute(const std::string&, TCallback, void*) {}

sqlite3* Database::open(const std::string&)
{
   return nullptr;
}

int Database::extract_callback(void*, int, char**, char**)
{
   return 1;
}

} // namespace database
} // namespace flexiblesusy

#endif
