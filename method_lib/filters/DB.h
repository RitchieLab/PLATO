// DB.h

#ifdef MYSQLDEFINED

#ifndef datahandling__DB_H__
#define datahandling__DB_H__

#include <mysql++.h>
#include "DataExcept.h"

using namespace std;

namespace Filters{

//typedef mysqlpp::Result DBResult;
typedef mysqlpp::Row DBRow;

///
/// Wrapper for API to connect to PLATO database in MySQL<br>
/// Database can be altered by wrapping a new API
///

/// Provides interaction with MySQL database
class DB{

  public:
    DB();
    DB(string hostname, string username, string pass, string databasenames);
    ~DB();
    
    void connect();
    void disconnect();  
    bool test_connection();
    void query(string query_string);
    bool execute(string mod_string);
    unsigned int get_last_id();
    
    DBRow fetch_row();
    //DBRow fetch_row(unsigned int index);
    string quote(string query_string);

    bool insert_ind(string dataset_id, unsigned int status, string all_1, 
        string all_2, string missing);
        
    void set_connection_params(string hostname, string username, string pass, unsigned int portnum, 
      string databasename);

    static bool EmptyRow;
       
  private:
    void initialize(string hostname, string username, string pass, unsigned int portnum, string databasename);
    
    mysqlpp::Connection conn;
    string host, user, passwd, dbname;
    unsigned int port;   
    mysqlpp::UseQueryResult dataset;
    
};

}
#endif

#endif
