//AlleleFrequencyFilter.cpp

#include "AlleleFrequencyFilter.h"

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
AlleleFrequencyFilter::AlleleFrequencyFilter():Filter("AlleleFreq"){
  initialize();
}


///
/// Alternative constructor
/// @param filterName name of filter
///
AlleleFrequencyFilter::AlleleFrequencyFilter(string filterName):Filter(filterName){
  initialize();
}

///
/// initializes maps used in switches for setting parameters
/// @return none
///
void AlleleFrequencyFilter::initialize(){
  freqThreshold = 0.0;

  ConfigMap["THRESHOLD"] = Threshold;

  FreqMap["ALL"] = "-all";
  FreqMap["ALLCHILDREN"] = "-all-children";
  FreqMap["RANDOMCHILDREN"] = "-random-child";
  FreqMap["UNAFFSPOUSESONLY"] = "-unaff-spouses";
  FreqMap["UNKSPOUSES"] = "-unk-spouses";
  
  frequency_option = "-all";
}



///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate AlleleFrequencyFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;

  start = clock();
  freq_calculator.resetDataSet(&dataset);
  end = clock();
  
  clock_t reset_time = end - start; 
 
  size_t num_run = resultList.size();

  start = clock();
  analyze(resultList, dataset);
  end = clock();

  time_estimate.seconds = double(end - start - reset_time) / (CLOCKS_PER_SEC * num_run);
  time_estimate.seconds *= num_models;
  
  // calculate total that will be passed through
  time_estimate.num_models = num_models;
  
  return time_estimate;
}



///
/// Runs marker genotype efficiency on dataset passed <br>
/// @param resultList ResultSet that may contain results from other filters 
/// and will contain results of this sequential replication analysis
/// @param dataset DataSet to be analyzed
/// @return none
/// @throws FilterExcept when no sequential replication type has been set
///
void AlleleFrequencyFilter::analyze(ResultSet & resultList, DataSet & dataset){ 

  // calculate position for p-value or likelihood ratio
  // will always be next value
  unsigned int threshold_index = 0;
  if(resultList.size() > 0){
    threshold_index = resultList.begin()->analysisScores.size();
  }
  bool erase_result = false;
  float result_value;

  freq_calculator.resetDataSet(&dataset);

  for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 
    erase_result = false;
    
    freq_calculator.calculate(currResult->genoCombination[0]);
    
    result_value = freq_calculator.getAone_freq()/
      (freq_calculator.getAtwo_freq()+freq_calculator.getAone_freq());
      
    // finding minor allele frequency
    if(result_value > 0.5)
      result_value = 1-result_value;
      
    if(result_value >= freqThreshold){
      currResult->analysisScores.push_back(result_value);
      ++currResult;
    }
    else
      resultList.erase(currResult++);
  } 

}



///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void AlleleFrequencyFilter::set_params(PARAMS params, DataSet* set){

  StepOptions options;
  
//   std::map<string, string>::iterator configIter;
  PARAMS::iterator configIter;
  string value;
  
  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){
    switch(ConfigMap[configIter->first]){
      case Threshold:
        freqThreshold = Stringmanip::stodouble(configIter->second); 
        if(freqThreshold < 0 || freqThreshold > 100){
          throw FilterExcept(configIter->first + 
            " must be greater than or equal to zero and less than or equal to 100 for " + 
            name + " filter\n\n");
        }
        break;
      case IncludeType: 
        value = Stringmanip::to_upper(configIter->second);
        if(FreqMap.find(value) != FreqMap.end()){
          frequency_option = FreqMap[value];
        }
        else{
          throw FilterExcept(configIter->second + " is not a valid option for paramater " +
            configIter->first + " for " + name + " filter\n\n");
        }
        break;    
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  }  
  
  options.setUp(frequency_option);
  freq_calculator.resetDataSet(set);
  freq_calculator.setOptions(options);
  
}

}




