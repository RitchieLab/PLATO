//BioFilterCombinations.cpp

#include "BioFilterCombinations.h"
#include <sstream>
#include "Stringmanip.h"


namespace Filters{
///
/// Constructor sets connection parameters
///
BioFilterCombinations::BioFilterCombinations(){
  
  string hostname = "munster.mc.vanderbilt.edu";
  string username = "plato";
  string pass = "genomewide";
  unsigned int portnum = 3306;
  string dbname = "Biofilter";
  
  // query connection
  combo_interval = 30000;
      
  conn.set_connection_params(hostname, username, pass, portnum, dbname);
  
  conn.connect(); 

}



BioFilterCombinations::BioFilterCombinations(string hostname, string username, string pass,
  string schema){
  
  unsigned int portnum = 3306;
   
   // query connection
  combo_interval = 30000;    
  conn.set_connection_params(hostname, username, pass, portnum, schema);
  
  conn.connect();

}


///
/// Destructor
///
BioFilterCombinations::~BioFilterCombinations(){
//  conn.disconnect();
}


///
/// Adds results to set based on limits 
/// 
void BioFilterCombinations::add_combos(){

  string sql = startsql + Stringmanip::itos(curr_res) + "," + increment;
  conn.query(sql);
  curr_res += combo_interval;

}


///
/// Queries database for list of snp combinations
/// @param tablename view or table to query
/// @param combination_size number of snps in combination
/// 
void BioFilterCombinations::fill_combinations(string tablename, unsigned int combination_size){


  // get total number of rows in table
  vector<string> snp_columns;
  
  stringstream ss;
  unsigned int i;
  
  for(i=1; i<= combination_size; i++){
    ss << "snp" << i;
    snp_columns.push_back(ss.str());
    ss.str("");
  }

  
  // current start result
  curr_res = 0;

  combo_size = combination_size;
  
  startsql = "SELECT concat('rs', cast(" + snp_columns[0] + " as char))";
  for(i=1; i<combination_size; i++){
    startsql += ", concat('rs', cast(" + snp_columns[i] + " as char))";
  }
  startsql += " FROM Biofilter." + tablename + " LIMIT ";
  
  increment = "30000";
}



///
/// Fills vector with as many combinations as are contained
/// in current view
/// @param combos two-dimensional vector containing the snp combinations
///
bool BioFilterCombinations::get_combinations(vector<vector<string> >& combos){
  bool query_finished = false;
  
  // add another increment of combos to results
  add_combos();
  
  vector<string> combo(combo_size,"");
  
  unsigned int i;

  DBRow r;
  
unsigned int total_fetched=0;  
  
  while(r = conn.fetch_row()){
    for(i=0; i<combo_size; i++){
      combo[i] = string(r[i]);
    }
    combos.push_back(combo);
total_fetched++;

  }
  
  if(combos.size() < combo_interval)
    query_finished = true;
  return query_finished;
}
}

