//ArmitageFilter.cpp
#include "ArmitageFilter.h"
#include <ChiSquare.h>

using namespace Methods;


namespace Filters{
///
/// constructor -- initialize variables
///
ArmitageFilter::ArmitageFilter():Filter("Armitage"){
  initialize();
}


///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
ArmitageFilter::ArmitageFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void ArmitageFilter::initialize(){
  use_raw_score = false;
  perm_finished = false;
  threshold = 1.0;
  n_perms = 1000;
  worst_score = 0.0;
  sort_priority = highest;
  total_type = ContingencyTable::Genotype;
  
  ConfigMap["THRESHOLD"] = Threshold;
  ConfigMap["USERAWSCORE"] = UseRawScore;
  ConfigMap["PERMUTATIONS"] = Permutations;
  ConfigMap["TEST"] = TotalTest;
  
  TotalTypeMap["ADDITIVE"] = ContingencyTable::Additive;
  TotalTypeMap["GENOTYPE"] = ContingencyTable::Genotype;  
  
}


///
/// Calculates Uncertainty Coefficient on all combinations in Result list
/// removes those that don't meet cut-off threshold
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void ArmitageFilter::analyze(ResultSet & resultList, DataSet & dataset){
  
//     if(!use_raw_score && !perm_finished) 
//       run_permutations(dataset);
    
    ContingencyTable table;
    table.set_current_totals(total_type);
    
    vector<Marker*>* markers = dataset.get_markers();
    int locus;
    
    for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 
//       table.get_counts(currResult->genoCombination[0], dataset);
      locus = (*markers)[currResult->genoCombination[0]]->getLoc();
      table.get_counts(locus, &dataset);
      table.set_current_totals(total_type);
      
      arm_score = call_armitage(table);
      // set score based on this table

      currResult->analysisScores.push_back(arm_score);
      
      if(currResult->analysisScores.back() < 0)
        currResult->analysisScores.back() = 0;

      // when not keeping raw scores -- use the permutation distribution 
      // to set the p value
      if(!use_raw_score){
        // determine degrees of freedom
        int num_rows = table.num_rows();
        if(num_rows != 1)
          --num_rows;
        int num_cols = table.num_cols();
        if(num_cols != 1)
          --num_cols;

        int df = num_cols * num_rows;
      
//        currResult->analysisScores.back() = perm.get_p_value(currResult->analysisScores.back());
        currResult->analysisScores.back() = ChiSquare::pfromchi(currResult->analysisScores.back(), df);
        if(currResult->analysisScores.back() > threshold){
          resultList.erase(currResult++);
        }
        else{
          ++currResult;
        }
      }
      else{
        if(currResult->analysisScores.back() < threshold){
          resultList.erase(currResult++);
        }
        else{
          ++currResult;
        }
      }   
  }  
}

///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate ArmitageFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;
 
  size_t num_run = resultList.size();
 
  start = clock();
  analyze(resultList, dataset);
  end = clock();

  time_estimate.seconds = double(end - start) / (CLOCKS_PER_SEC * num_run);
  time_estimate.seconds *= num_models;
  
  // calculate total that will be passed through
  // depends on threshold (converted to p value if not already)
  float pthresh = threshold;
  if(use_raw_score){
    int df = 2; // usual degrees of freedom based on 3 X 2 table
    ChiSquare::pfromchi(threshold, df);
  }
  time_estimate.num_models = num_models * pthresh;
  
  return time_estimate;
}


///
/// Converts totals from ContingencyTable into appropriate
/// vector for running armitage test and runs test
/// @param table ContingencyTable
///
float ArmitageFilter::call_armitage(ContingencyTable& table){

  vector<vector<float> > float_totals = table.get_vector();
  
  vector<int> temp(float_totals[0].size(),0);
  vector<vector<int> > int_totals(float_totals.size(), temp);
  
  unsigned int i,j;
  
  for(i=0; i<float_totals.size(); i++){
    for(j=0; j<float_totals[i].size(); j++){
      int_totals[i][j] = int(float_totals[i][j]);
    }
  }
  
  return arm_calculator.armitage(int_totals);
}


///
/// Calculates likelihood ratio on table
/// @param res Result containing genotype to check
/// @param table ContingencyTable containing totals to use in calcuation
/// 
void ArmitageFilter::analyze(Result& res, ContingencyTable& original_table){
  
//   or_calculator.calculate(&original_table);
  original_table.set_current_totals(total_type);  
  res.analysisScores.push_back(call_armitage(original_table));
//   res.analysisScores.push_back(or_calculator.getOddsRatio());
}


///
/// Uses Permutations object to run permutations on the dataset
/// @param dataset DataSet 
/// @return
///
void ArmitageFilter::run_permutations(DataSet& dataset){
  bool original_use_raw_score = use_raw_score;
  use_raw_score = true;
  float orig_threshold = threshold;
  threshold = 0.0;
  perm.add_filter(this);
  perm.run_permutations(dataset, n_perms);
  use_raw_score = original_use_raw_score;
  threshold = orig_threshold;
}


///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void ArmitageFilter::set_params(PARAMS params, DataSet* set){
  std::map<string, string>::iterator configIter;
  string value;
 
  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){ 
    switch(ConfigMap[configIter->first]){
      case UseRawScore:
        value = Stringmanip::to_upper(configIter->second);
        if(value.find("TRUE") != string::npos)
          use_raw_score = true;
        else if(value.find("FALSE") != string::npos)
          use_raw_score = false;
        else
          throw FilterExcept(configIter->first + " must be either True or False for " +
            name + " filter\n\n");
        break;
      case Threshold:
        threshold = Stringmanip::stodouble(configIter->second);
        if(threshold < 0){
          throw FilterExcept(configIter->first + " must be greater than zero for " + 
            name + " filter\n\n");          
        }
        break;
      case TotalTest:
        value = Stringmanip::to_upper(configIter->second);
        if(TotalTypeMap.find(value) != TotalTypeMap.end()){
          total_type = TotalTypeMap[value];
        }
        else{
          throw FilterExcept(configIter->second + " is not a valid option for paramater " +
            configIter->first + " for " + name + " filter\n\n");
        }        
        break;  
      case Permutations:
        n_perms = Stringmanip::stouint(configIter->second);
        if(n_perms < 0 || n_perms > 100000)
          throw FilterExcept(configIter->first + " must be between 0 and 100000 for " +
            name + " filter\n\n");
        break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  }  
 
  switch(total_type){
    case ContingencyTable::Genotype:
      name += ":Genotype";
      break;     
    case ContingencyTable::Additive:
      name += ":Additive";
      break;
    case ContingencyTable::Allele:
    case ContingencyTable::Dominant:
    case ContingencyTable::Recessive:
      break;
  };
  
}
}
