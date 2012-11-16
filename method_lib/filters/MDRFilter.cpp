//MDRFilter.cpp

#include "MDRFilter.h"

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
MDRFilter::MDRFilter():Filter("MDR Filter"){
  name = "MDR";
  initialize();
}

///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
MDRFilter::MDRFilter(string filterName):Filter(filterName){
  initialize();
}

///
/// initializes maps used in switches for setting parameters
/// @return none
///
void MDRFilter::initialize(){
//   sort_priority = highest;
//   worst_score = 0.0;
//   n_perms = 0;
//   perm_completed = false;

  accThreshold = 0.50;
  ConfigMap["MODTHRESHCALC"] = ModelThreshCalculated;
  ConfigMap["ACCTHRESHOLD"] = Threshold;
  
//   ConfigMap["PERMUTATIONS"] = Permutations;
  
}



///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @param 
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void MDRFilter::set_params(PARAMS params, DataSet* dataset){
  std::map<string, string>::iterator configIter;
  string value;
  
  StepOptions options;
  
  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){
    switch(ConfigMap[configIter->first]){
      case Threshold:
        accThreshold = Stringmanip::stodouble(configIter->second);
        if(accThreshold < 0){
          throw FilterExcept(configIter->first + " must be greater than zero for " + 
            name + " filter\n\n");          
        }
        break;
      case ModelThreshCalculated:
        value = Stringmanip::to_upper(configIter->second);
        if(value.compare("SET") == 0 || value.compare("ON") == 0){
          options.setOnlySetThreshold(true);
        }
        else if(value.compare("MODEL")==0){
          options.setOnlySetThreshold(false);
        }
        else{
          throw FilterExcept(configIter->second + " is not a valid option for parameter " +
            configIter->first + " for " + name + " filter\n\n");
        }
        break;
//       case Permutations:
//         n_perms = Stringmanip::stouint(configIter->second);
//         if(n_perms < 1 || n_perms > 100000)
//           throw FilterExcept(configIter->first + " must be between 1 and 100000 for " +
//             name + " filter\n\n");
//         break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  } 
  mdr_calculator.set_parameters(&options);
}



///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate MDRFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;
  
  start = clock();
  mdr_calculator.resetDataSet(&dataset);
  end = clock();

  clock_t reset_time = end - start;
 
  size_t num_run = resultList.size();
 
  start = clock();
  analyze(resultList, dataset);
  end = clock();

  time_estimate.seconds = double(end - start - reset_time) / (CLOCKS_PER_SEC * num_run);
  time_estimate.seconds *= num_models;
  
  // calculate total that will be passed through
  time_estimate.num_models = (1-get_model_fraction(model_size)) * num_models;
  
  return time_estimate;
}


///
/// Returns fraction of models that will meet threshold based on 
/// model size.  Distributions were created with null data.
/// @param model_size int
/// @return fraction that will meet or exceed threshold
///
double MDRFilter::get_model_fraction(int model_size){
  
  vector<double> standard_dev;
  standard_dev.push_back(0);
  standard_dev.push_back(0.00640);
  standard_dev.push_back(0.01027);
  standard_dev.push_back(0.01167);
  standard_dev.push_back(0.01514);
  
  vector<double> mean;
  mean.push_back(0);  // means appear to increase by about .5 of the difference from last
  mean.push_back(0.51486);
  mean.push_back(0.52883);
  mean.push_back(0.54859);
  mean.push_back(0.57855);
  // predict 5 way models should have a mean of around 0.6235
  
  // if greater than max use max size as estimate
  if(model_size > 4){
    model_size = 4;
  }
  
  return normdist(accThreshold, mean[model_size], standard_dev[model_size]);

}



///
/// Runs MDR on dataset passed <br>
/// @param resultList ResultSet that may contain results from other filters 
/// and will contain results of this MDR analysis
/// @param dataset DataSet to be analyzed
/// @return none
///
void MDRFilter::analyze(ResultSet & resultList, DataSet & dataset){ 

  mdr_calculator.resetDataSet(&dataset);

  for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 

     mdr_calculator.calculate(currResult->genoCombination);
    
     currResult->analysisScores.push_back(mdr_calculator.getBalancedAcc());
     
     if(currResult->analysisScores.back() < 0)
      currResult->analysisScores.back() = 0;
     
    if(currResult->analysisScores.back() < accThreshold){
      resultList.erase(currResult++);
    }
    else 
      currResult++;
    
  } 

}
}

