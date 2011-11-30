//LikelihoodRatioFilter.cpp
#include "LikelihoodRatioFilter.h"

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
LikelihoodRatioFilter::LikelihoodRatioFilter():Filter("LikelihoodRatio"){
  initialize();
}


///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
LikelihoodRatioFilter::LikelihoodRatioFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void LikelihoodRatioFilter::initialize(){
  use_raw_score = false;
  perm_finished = false;
  threshold = 1.0;
  n_perms = 1000;
  worst_score = 0.0;
  sort_priority = highest;
  
//   total_type = ContingencyTable::Allele;
  
  totalType = "ALLELE";
  
  ConfigMap["THRESHOLD"] = Threshold;
  ConfigMap["USERAWSCORE"] = UseRawScore;
  ConfigMap["PERMUTATIONS"] = Permutations;
  ConfigMap["TEST"] = TotalTest;
  
  
  NullMean["ADDITIVE"]=1.98205565593678;
  NullMean["ALLELE"]=0.996689821221706;
  NullMean["GENOTYPE"]=1.98205565593678;
  NullMean["DOMINANT"]=0.993146095746896;
  NullMean["RECESSIVE"]=0.987662397101939;
  NullStdDev["ADDITIVE"]=1.9678949619182;
  NullStdDev["ALLELE"]=1.4088947820973;
  NullStdDev["GENOTYPE"]=1.9678949619182;
  NullStdDev["DOMINANT"]=1.40510479800273;
  NullStdDev["RECESSIVE"]=1.3840415409424;
 
}


///
/// Calculates Uncertainty Coefficient on all combinations in Result list
/// removes those that don't meet cut-off threshold
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void LikelihoodRatioFilter::analyze(ResultSet & resultList, DataSet & dataset){
    llr_calculator.resetDataSet(&dataset);
  
    if(!use_raw_score && !perm_finished) 
      run_permutations(dataset);  
    for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 
//       table.get_counts(currResult->genoCombination[0], dataset);
      llr_calculator.calculate(currResult->genoCombination[0]);
      // set score based on this table
//       analyze(*currResult, table);

      currResult->analysisScores.push_back(llr_calculator.getLikelihoodRatio());
      
      if(currResult->analysisScores.back() < 0)
        currResult->analysisScores.back() = 0;     
      
      // when not keeping raw scores -- use the permutation distribution 
      // to set the p value
      if(!use_raw_score){
        currResult->analysisScores.back() = perm.get_p_value(currResult->analysisScores.back());
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
Filter::ProcessEstimate LikelihoodRatioFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;

  start = clock();
  llr_calculator.resetDataSet(&dataset);
  end = clock();
  
  clock_t reset_time = end - start;
 
  size_t num_run = resultList.size();
 
  unsigned int original_perms = n_perms;
  n_perms = 0;
  start = clock();
  analyze(resultList, dataset);
  end = clock();
  n_perms = original_perms;
  time_estimate.seconds = double(end - start - reset_time) / (CLOCKS_PER_SEC * num_run);
  time_estimate.seconds *= num_models * (n_perms+1);
  
  // calculate total that will be passed through
  if(n_perms <= 0 || use_raw_score)
    time_estimate.num_models = (1-get_model_fraction(model_size)) * num_models;
  else
    time_estimate.num_models = threshold * num_models;
  
  return time_estimate;
}


///
/// Returns fraction of models that will not meet threshold based on 
/// model size.  Distributions were created with null data.
/// @param model_size int
/// @return fraction that will meet or exceed threshold
///
double LikelihoodRatioFilter::get_model_fraction(int model_size){
  // return all when threshold is zero
  if(threshold <= 0){
    return 0.0;
  }
  double standard_dev = NullStdDev[totalType];
  double mean = NullMean[totalType];
  
  return normdist(threshold, mean, standard_dev);
}

///
/// Calculates likelihood ratio on table
/// @param res Result containing genotype to check
/// @param table ContingencyTable containing totals to use in calcuation
/// 
void LikelihoodRatioFilter::analyze(Result& res, ContingencyTable& original_table){
  
  llr_calculator.calculate(&original_table);
  res.analysisScores.push_back(llr_calculator.getLikelihoodRatio());
  
}


///
/// Uses Permutations object to run permutations on the dataset
/// @param dataset DataSet 
/// @return
///
void LikelihoodRatioFilter::run_permutations(DataSet& dataset){
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
void LikelihoodRatioFilter::set_params(PARAMS params, DataSet* set){
  std::map<string, string>::iterator configIter;
  string value;
  StepOptions options;
 
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
        totalType = Stringmanip::to_upper(configIter->second);
        options.setLikelihoodRatioTotalType(totalType);
//         if(TotalTypeMap.find(value) != TotalTypeMap.end()){
//           total_type = TotalTypeMap[value];
//         }
//         else{
//           throw FilterExcept(configIter->second + " is not a valid option for paramater " +
//             configIter->first + " for " + name + " filter\n\n");
//         }        
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
 
  name += ":" + Stringmanip::to_upper(totalType);
  
  llr_calculator.set_parameters(&options);
  
}
}
