//BioRecordResults.h

#ifndef __BIORECORDRESULTS_H__
#define __BIORECORDRESULTS_H__

#include "DB.h"
#include "ResultSet.h"
#include <DataSet.h>

namespace Filters{
/// Output results to database
class BioRecordResults{

  public:
    /// Contructor
    BioRecordResults();
    
    /// Destructor
    ~BioRecordResults();
       
    /// Record results in designated table
    void record_results(ResultSet & resultset, 
      Methods::DataSet & dataset, string table_name,  bool interact_included, bool include_cov);
    
    /// Set connection parameters for this recorder of results
    void set_connect_params(string hostname, string username, string pass,
      string dbname, uint portnum=3306);
    
    /// Establishes connection 
    void connect();
    
    /// Creates table if doesn't exist for output
    void check_table(int maxModelSize, bool interact_included, string table_name, string schema_name);
    
  private:
    
   /// Get rs numbers for the indicated loci    
   void get_rs_numbers(vector<string>& rs_nums, vector<uint> & indexes,
    Methods::DataSet& dataset);   

   int num_coefficients(unsigned int size, bool interact_included);

    DB conn;

};

}

#endif
