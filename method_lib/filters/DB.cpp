// DB.cpp

#ifdef MYSQLDEFINED

#include "DB.h"
namespace Filters{
bool DB::EmptyRow = false;

///
///Constructor establishes connection to database
///
DB::DB(){

  host="munster.mc.vanderbilt.edu"; // mac pro
  user = "plato";
  passwd = "genomewide"; 
  port = 3306;
  dbname = "Analysis";

  initialize(host, user, passwd, port, dbname);
//  connect();
}


///
/// Alternate constructor allowing for connection to alternative database
/// @param hostname
/// @param username
/// @param pass
/// @param databasename
///
DB::DB(string hostname, string username, string pass, string databasename){
  initialize(hostname, username, passwd, 3306, databasename);
}


///
/// Tests the connection 
///
bool DB::test_connection(){
  return conn.connected();
}


///
///Destructor terminates connection to database
///
DB::~DB(){
  if(conn.connected())
    disconnect();
}


///
///Connects with database
///@return 
///
void DB::connect(){
  conn.connect(dbname.c_str(), host.c_str(), user.c_str(), passwd.c_str(), port);
}


///
///Closes connection 
///@return 
///
void DB::disconnect(){
  conn.disconnect();
}


///
///Runs query passed to it <br>
///The results can then be extracted per row
///using the fetch_row method
///@param query_string string consisting of sql query
///@return none
///@throws exception when fails
///
void DB::query(string query_string){
  mysqlpp::Query current_query = conn.query();
  current_query << query_string;
  try{
    dataset = current_query.use();
  }
  catch (const mysqlpp::BadQuery& er) {
    // Handle any query errors
    string errormess(er.what()); 
    throw DataExcept("Query error: " + errormess);
  }
  catch(const mysqlpp::Exception& er){
    // Handle any other errors
    throw DataExcept(query_string + " failed to execute: " + er.what());
  }
}


///
/// Used when modifying data <br>
/// If fails throws exception.
/// @param mod_string string containing SQL for 
///        update or insertion
/// @return true for success and false for failure
/// @throws DataExcept for failure
///
bool DB::execute(string mod_string){
  mysqlpp::Query current_query = conn.query();
  bool query_success;
  
  try{
    query_success = current_query.exec(mod_string);
  }
  catch (const mysqlpp::BadQuery& er) {
    // Handle any query errors
    string errormess(er.what());
    throw DataExcept("Query error: " + errormess);
  }
  catch(const mysqlpp::Exception& er){
    // Handle any other errors
    throw DataExcept(mod_string + " failed to execute: " + er.what());
  }
  return query_success;
}


///
/// Returns last auto_incremented ID for session <br>
/// If fails throws exception.
/// @return ID 
/// @throws DataExcept on error
///
unsigned int DB::get_last_id(){
  string sql = "SELECT Last_insert_id()";
  mysqlpp::Query current_query = conn.query();
  current_query << sql;
  try{
    dataset = current_query.use();
  }
  catch(const mysqlpp::Exception& er){
    // Handle any other errors
    throw DataExcept(sql + " failed to execute: " + er.what());
  }
  
  DBRow r = dataset.fetch_row();
  unsigned int id = r[int(0)];
//  dataset.purge();
  return id;
}


///
///Quotes the string for use in a query
///and escapes any quotes within the string
///@param query_string String containing text to quote for query
///@return new quoted string
///
string DB::quote(string query_string){
  stringstream ss;
  ss << mysqlpp::quote << query_string;
  return ss.str();
}


///
///Fetches next row in the dataset<br>
///When no row left to fetch an empty row 
///is returned
///@return next row in the set
///
DBRow DB::fetch_row(){
  
  try{
    return dataset.fetch_row();
  }
  catch(const mysqlpp::Exception& er){
    DBRow empty_row;
    return empty_row;
  }
}


///
/// Inserts individual information as 3 Blobs in the MySQL database
/// @param dataset_id string listing database ID for this individual
/// @param status Individual affection status (0 or 1)
/// @param all_1 bit string of alleles
/// @param all_2 bit string of alleles
/// @param missing bit string with bit on whenever data is missing
/// @return true if successful
///
bool DB::insert_ind(string dataset_id, unsigned int status, string all_1, 
  string all_2, string missing){
  ostringstream strbuf;  
  try{
    mysqlpp::Query query = conn.query();
    
    strbuf << "INSERT INTO inds(dataset_id, status, all_1, all_2, missing) VALUES(" << 
      dataset_id << "," << status << ",\"" <<
      mysqlpp::escape << all_1 << "\"," << mysqlpp::escape << all_2
      << "\"," << mysqlpp::escape << missing << "\")";
      query.exec(strbuf.str());
  }
    catch(const mysqlpp::Exception& er){
    // Handle any other errors
    throw DataExcept(strbuf.str() + " failed to execute: " + er.what());
  }
  return true;    
}

///
/// Sets connection parameters
/// @param hostname
/// @param username
/// @param pass
/// @param portnum
///
void DB::set_connection_params(string hostname, string username, string pass, unsigned int portnum, 
  string databasename){
  
  host=hostname; // mac pro
  user = username;
  passwd = pass; 
  port = portnum;
  dbname = databasename;  
  
}


///
///Sets initial values for connection parameters
///@return 
///
void DB::initialize(string hostname, string username, string pass, unsigned int portnum, string databasename){

  set_connection_params(hostname, username, pass, portnum, databasename);
    
}
}

#endif
