//BioRecordResults.cpp

#include "BioRecordResults.h"
#include "Stringmanip.h"
#include <math.h>
#include <ComboGenerator.h>

using namespace Methods;
namespace Filters{
///
/// Constructor sets connection parameters
///
BioRecordResults::BioRecordResults(){
  string hostname = "munster.mc.vanderbilt.edu";
  string username = "plato";
  string pass = "genomewide";
//   unsigned int portnum = 3306;
  string dbname = "Analysis";
}

///
/// Destructor
///
BioRecordResults::~BioRecordResults(){
  if(conn.test_connection())
    conn.disconnect();
}


///
/// Sets up connection
/// 
void BioRecordResults::connect(){
  conn.connect();
}


///
/// Set connection parameters
/// @param hostname
/// @param username
/// @param pass
/// @param dbname
/// @param portnum
///
void BioRecordResults::set_connect_params(string hostname, string username, string pass,
  string dbname, unsigned int portnum){
  conn.set_connection_params(hostname, username, pass, portnum, dbname); 
}


///
/// Records all results in the specified set.
/// @param resultset ResultSet with results to insert
/// @param dataset DataSet with loci information
/// @param table_name name of database table results will be written to
/// @param interact_included when true the interaction term was included
/// @param include_cov when true covariates are included to be inserted into database
/// @return 
/// @throws PlatoExcept on error
///
void BioRecordResults::record_results(ResultSet& resultset,
  DataSet& dataset, string table_name, bool interact_included, bool include_cov){
  
  unsigned int data_fields = resultset.begin()->analysisScores.size();

  // need to bunch results to have all the same size ones grouped together
  // so insert as batches of results
  ResultIter currRes = resultset.begin();
  
  string sql;
  
  unsigned int resultsInQuery=0, resultsMax=10000;
  
  while(currRes != resultset.end()){
  
    // build starting string based on first result to insert
    unsigned int nLoci = currRes->genoCombination.size();
    
    string startsql = "INSERT INTO " + table_name + " (";
    
    stringstream ss;
    ss << "snp" << 1;
    for(unsigned int snp=2; snp <= nLoci; snp++){
      ss << ",snp" << snp;
    }
    
    startsql += ss.str();
    startsql += ",likelihood_ratio";
    
    // add coefficients
    ss.str("");
    int numcoeff = num_coefficients(nLoci, interact_included);
    for(int i=1; i<numcoeff; i++){
      ss << ",coeff" << i;
    }
    ss << ",coeff" << numcoeff;
    startsql += ss.str() + ")VALUES";

    sql = startsql;
    
    vector<string> rs_numbers;
    
    // need to build on to this result with all results 
    // that are the same size until reach the end of the list
    // or reach the max size
    while(currRes->genoCombination.size() == nLoci && currRes != resultset.end() &&
      resultsInQuery <= resultsMax){
    
      get_rs_numbers(rs_numbers, currRes->genoCombination, dataset);
      
      if(resultsInQuery>0)
        sql +=",";
      
      sql += "(";
      
      for(unsigned int i=0; i<rs_numbers.size(); i++){
        sql += rs_numbers[i] + ",";
      }
      
      stringstream ss;
      
      // add likelihood score and coefficients
      ss << currRes->analysisScores[0] << ",";
      for(int i=1; i<numcoeff; i++){
        ss << currRes->analysisScores[i] << ",";
      }
      ss << currRes->analysisScores[numcoeff] << ")";
      
      sql += ss.str();
      
    
      resultsInQuery++;
      currRes++;
    }  
    
    // if any results in query submit
    if(resultsInQuery > 0){
      while(!conn.test_connection())
        conn.connect();
        conn.execute(sql); 
    }
    
    
    resultsInQuery = 0;
  }
  
}


///
/// Creates table if doesn't exist for output
/// @param maxModelSize 
/// @param interact_included
/// @param table_name
/// @param schema_name
///
void BioRecordResults::check_table(int maxModelSize, bool interact_included, string table_name, string schema_name){
  
  string sql = "select count(*) from information_schema.tables where table_schema = '" + schema_name 
    + "' and table_name = '" + table_name + "'";
    
  while(!conn.test_connection())
    conn.connect();

  
  conn.query(sql);
  
  DBRow r = conn.fetch_row();
  
  if(int(r[0]) == 1){
    conn.disconnect();
    return;
  }
  
  conn.disconnect();
  
  conn.connect();
  // if reached here need to create the table
  
  sql = "CREATE TABLE " + schema_name + "." + table_name + "(";
  
  stringstream ss;
  for(unsigned int i=1; i<= maxModelSize; i++){
    ss << "snp" << i << " INT,";
  }
  
  // add likelihood ratio
  ss << "likelihood_ratio FLOAT,";
  
  // need to know how many coefficients to include
  int numCoeff = num_coefficients(maxModelSize, interact_included);
  
  for(unsigned int i=1; i<numCoeff; i++){
    ss << "coeff" << i << " FLOAT,";
  }
  ss << "coeff" << numCoeff << " FLOAT)";
  
  sql += ss.str();

  conn.execute(sql);
  
}


///
/// returns expected number of coefficients based on model size and whether
/// there are interactions included
/// @param res Result
/// @param interact_included 
/// @return number of coefficients
///
int BioRecordResults::num_coefficients(unsigned int size, bool interact_included){
  
  unsigned int nLoci = size;
  
  if(!interact_included)
    return nLoci;
  
  int num_coeff = nLoci;
  
  if(num_coeff > 1){
    ComboGenerator combo;
  
    for(unsigned int i=2; i <= nLoci ; i++){
      num_coeff += int(combo.calc_combinations(nLoci, i));
    }
  }
  
  return num_coeff;
  
}


///    
/// Records all results in the specified set in the database.<br>
/// Order of analyis must be set using set_analysis function.
/// @param resultset ResultSet with results to insert
/// @param dataset DataSet with loci information
/// @param table_name name of database table results will be written to
/// @param interact_included when true the interaction term was included
/// @param include_cov when true covariates are included to be inserted into database
/// @return 
/// @throws PlatoExcept on error
///
// void BioRecordResults::record_results(ResultSet & resultset, 
//   DataSet & dataset, string table_name, bool interact_included, bool include_cov){
// 
//     ////////////////////////////////////
//     // covariate matrix is as follows  
//     //      snp1   snp2    Int
//     // B1   cov1
//     // B2   cov2   cov3
//     // B3   cov4   cov5    cov6
//     /////////////////////////////////////
//   // use amount of data stored 
//   unsigned int data_fields = resultset.begin()->analysisScores.size();
// // cout << "number of data fields = " << data_fields << endl;
//   unsigned int final_field_index = data_fields - 1;
//   
//   string startsql = "INSERT INTO " + table_name + " (snp1, snp2, likelihood_ratio, coeff1, coeff2";
//   if(interact_included){
//     startsql += ", coeff3";
//   }
//   
//   if(include_cov){
//     startsql += ", cov1, cov2, cov3";
//     if(interact_included){
//       startsql += ", cov4, cov5, cov6";
//     }
//   }
//   
//   startsql += ")VALUES";
//   
//   unsigned int currValue, totalValues, resultsInQuery=0, resultsMax=10000, totalInserted;
//   string sql = startsql;
//   vector<string> rs_numbers;
//   string snp1, snp2;
// 
//   ResultIter resIter = resultset.begin();
//   
//   totalValues = resIter->analysisScores.size();
// 
//   sql = startsql;
//   
//   for(; resIter != resultset.end(); resIter++){
//     get_rs_numbers(rs_numbers, resIter->genoCombination, dataset);
//     
//     snp1 = rs_numbers[0];
//     snp2 = rs_numbers[1];
//     
//     totalInserted = 0; 
//     if(resultsInQuery > 0)
// 	    sql += ",";    
// 
// 	  sql+= "(";  
// 	  
// 	  sql += snp1;
// 	  sql += "," + snp2 + ",";
// 	
//     // find first result that is valid for insertion
//     for(currValue=0; currValue < final_field_index; currValue++){
//         if(isinf(resIter->analysisScores[currValue]) || isnan(resIter->analysisScores[currValue])){
//           sql += "NULL,";
//         }
//         else
//           sql += Stringmanip::itos(resIter->analysisScores[currValue]) + ",";
//     }
//     if(isinf(resIter->analysisScores[currValue]) || isnan(resIter->analysisScores[currValue])){
//       sql += "NULL";
//     }
//     else
//       sql +=  Stringmanip::itos(resIter->analysisScores[currValue]);
//     
//     
//     sql += ")";
// 
//     ++resultsInQuery;
// 
//     if(resultsInQuery >= resultsMax){
//       while(!conn.test_connection())
//         conn.connect();
// //  cout << "executing sql " << endl << sql << endl;
// // exit(1);
//       conn.execute(sql);
// 
//       resultsInQuery=0;
//       sql = startsql;
//     }
//   }
//   if(resultsInQuery>0){
//     while(!conn.test_connection())
//       conn.connect();
// //  cout << "executing sql " << endl << sql << endl;
// // exit(1);
//     conn.execute(sql);
//   }
// // exit(1);
// }


///
/// Converts locus indexes into rs numbers
/// strips the rs if part of the string
/// @param rs_nums store the rs numbers as strings for insertion into database
/// @param indexes
/// @param dataset DataSet
///
void BioRecordResults::get_rs_numbers(vector<string>& rs_nums, vector<unsigned int> & indexes,
  DataSet& dataset){
  
  rs_nums.assign(indexes.size(), "");
  string string_id;
  string::size_type pos;
  
  vector<Marker*> markers = *(dataset.get_markers());
  vector<int> marker_map = *(dataset.get_marker_map());
  
  for(unsigned int i=0; i<indexes.size(); i++){
   //  string_id = dataset.get_locus_name(indexes[i]);
    
    string_id = markers[marker_map[indexes[i]]]->getRSID();
    
    pos = string_id.find("rs");
    // remove the 'rs' when found
    if(pos != string::npos){
      string_id.replace(pos, 2, "");
    }
    rs_nums[i] = string_id;
  
  }
}
}
