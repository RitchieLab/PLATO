//ChiSquareFilter.cpp
#include "ChiSquareFilter.h"

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
ChiSquareFilter::ChiSquareFilter():Filter("Chi Square"){
  initialize();
}


///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
ChiSquareFilter::ChiSquareFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void ChiSquareFilter::initialize(){
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
  
  TotalTypeMap["ALLELE"] = ContingencyTable::Allele;
  TotalTypeMap["GENOTYPE"] = ContingencyTable::Genotype;  
  
}


///
/// Calculates Uncertainty Coefficient on all combinations in Result list
/// removes those that don't meet cut-off threshold
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void ChiSquareFilter::analyze(ResultSet & resultList, DataSet & dataset){

    chi_calculator.resetData(&dataset);
    
    for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 
      currResult->analysisScores.push_back(call_chisquare(currResult->genoCombination[0]));


      // when not keeping raw scores -- use the permutation distribution 
      // to set the p value
      if(!use_raw_score){
        switch(total_type){
          case ContingencyTable::Genotype:
            currResult->analysisScores.back() = chi_calculator.getGenotypicPval();
          break;
          case ContingencyTable::Allele:
            currResult->analysisScores.back() = chi_calculator.getAllelicPval();
         case ContingencyTable::Additive:
         case ContingencyTable::Dominant:
         case ContingencyTable::Recessive:
          break;
        };
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
/// Converts totals from ContingencyTable into appropriate
/// vector for running armitage test and runs test
/// @param locus
///
float ChiSquareFilter::call_chisquare(int locus){

 chi_calculator.calculate(locus);
 
 float return_value = 1.0;
 
 switch(total_type){
    case ContingencyTable::Genotype:
      return_value = chi_calculator.getGenotypicChi();
      break;     
    case ContingencyTable::Allele:
      return_value = chi_calculator.getAllelicChi();
      break;
    case ContingencyTable::Additive:
    case ContingencyTable::Dominant:
    case ContingencyTable::Recessive:
    break;
  };
  return return_value;
}


///
/// Calculates chi square on table
/// @param res Result containing genotype to check
/// @param table ContingencyTable containing totals to use in calcuation
/// 
void ChiSquareFilter::analyze(Result& res, ContingencyTable& original_table){
  
  res.analysisScores.push_back(call_chisquare(res.genoCombination[0]));
}


///
/// Uses Permutations object to run permutations on the dataset
/// @param dataset DataSet 
/// @return
///
void ChiSquareFilter::run_permutations(DataSet& dataset){
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
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate ChiSquareFilter::estimate_run_time(double num_models, DataSet & dataset, 
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
    switch(total_type){
      case ContingencyTable::Genotype:
        df = 2;
        break;     
      case ContingencyTable::Allele:
        df = 1;
        break;
    case ContingencyTable::Additive:
    case ContingencyTable::Dominant:
    case ContingencyTable::Recessive:
    break;
    };
    ChiSquare::pfromchi(threshold, df);
  }
  time_estimate.num_models = num_models * pthresh;
  
  return time_estimate;
}



///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void ChiSquareFilter::set_params(PARAMS params, DataSet* set){
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
    case ContingencyTable::Allele:
      name += ":Allelic";
      break;
    case ContingencyTable::Additive:
    case ContingencyTable::Dominant:
    case ContingencyTable::Recessive:
      break;
  };
  
  // set dataset for running chi square
  chi_calculator.resetData(set);
  
}
}
