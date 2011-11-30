//MDRPDTFilter.cpp

#include "MDRPDTFilter.h"


using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
MDRPDTFilter::MDRPDTFilter():Filter("MDRPDT Filter"){
  name = "MDRPDT";
  initialize();
}

///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
MDRPDTFilter::MDRPDTFilter(string filterName):Filter(filterName){
  initialize();
}

///
/// initializes maps used in switches for setting parameters
/// @return none
///
void MDRPDTFilter::initialize(){

  threshold = 0.0;
  n_perms = 0;
  
  ConfigMap["THRESHOLD"] = Threshold;
  ConfigMap["RANDSEED"] = RandSeed;
  ConfigMap["PTESTS"] = Ptests;
  ConfigMap["CROSSVALS"] = Crossvals;
  ConfigMap["MINCOMBO"] = MinComb;
  ConfigMap["MAXCOMBO"] = MaxComb;
  
}



///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @param 
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void MDRPDTFilter::set_params(PARAMS params, DataSet* dataset){
  std::map<string, string>::iterator configIter;
  string value;
  
  StepOptions options;
  
  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){
    switch(ConfigMap[configIter->first]){
      case Threshold:
        threshold = Stringmanip::stodouble(configIter->second);
        if(threshold < 0){
          throw FilterExcept(configIter->first + " must be greater than zero for " + 
            name + " filter\n\n");          
        }
        break;
      case RandSeed:
        {
        int rand_seed = Stringmanip::stouint(configIter->second);
        if(rand_seed < 0)
          throw FilterExcept(configIter->first + " must be greater than 0 for " +
            name + " filter\n\n");
        options.setMDRPDTRandSeed(rand_seed);
        }
        break;
      case Ptests:
        {
        n_perms = Stringmanip::stouint(configIter->second);
        if(n_perms < 0 || n_perms > 100000)
          throw FilterExcept(configIter->first + " must be between 0 and 100000 for " +
            name + " filter\n\n");
        options.setMDRPDTNumPTests(n_perms);
        }
        break;
      case Crossvals:
        {
          int cv = Stringmanip::stouint(configIter->second);
          if(cv < 0 || cv > 100000)
            throw FilterExcept(configIter->first + " must be between 0 and 100000 for " +
              name + " filter\n\n");
          options.setMDRPDTNumCrossVals(cv);          
        }
        break;
      case MinComb:
        {
        int mincomb = Stringmanip::stouint(configIter->second);
          if(mincomb < 1 || mincomb > 10)
            throw FilterExcept(configIter->first + " must be between 1 and 10 for " +
              name + " filter\n\n");
          options.setMDRPDTMinCombo(mincomb);   
        }
        break;
      case MaxComb:
        {
        int maxcomb = Stringmanip::stouint(configIter->second);
          if(maxcomb < 1 || maxcomb > 10)
            throw FilterExcept(configIter->first + " must be between 1 and 10 for " +
              name + " filter\n\n");
          options.setMDRPDTMaxCombo(maxcomb);   
        }
        break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  } 
  mdrpdt_calculator.set_parameters(&options);
}


///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate MDRPDTFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;
  
  StepOptions options;
  unsigned int orig_ptests = options.getMDRPDTNumPTests();
  n_perms = 0;
  options.setMDRPDTNumPTests(n_perms);

  start = clock();
  mdrpdt_calculator.resetDataSet(&dataset);

  end = clock();

  clock_t reset_time = end - start;

 
  size_t num_run = resultList.size();

  start = clock();
  analyze(resultList, dataset);
  end = clock();

  time_estimate.seconds = double(end - start - reset_time) / (CLOCKS_PER_SEC * num_run *(1+orig_ptests));
  time_estimate.seconds *= num_models;
  
  // calculate total that will be passed through
  if(orig_ptests > 0){
    time_estimate.num_models = threshold * num_models;
  }
  else
    time_estimate.num_models = num_models;

  
  options.setMDRPDTNumPTests(orig_ptests);
  n_perms = orig_ptests;
  
  return time_estimate;
}


///
/// Runs MDR on dataset passed <br>
/// @param resultList ResultSet that may contain results from other filters 
/// and will contain results of this MDR analysis
/// @param dataset DataSet to be analyzed
/// @return none
///
void MDRPDTFilter::analyze(ResultSet & resultList, DataSet & dataset){ 

  mdrpdt_calculator.resetDataSet(&dataset);

  vector<float> t_stats, t_testing, mor_vals;
  vector<float>::iterator test_iter, t_iter, m_iter;

  for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){ 

     mdrpdt_calculator.calculate(currResult->genoCombination);

     t_stats = mdrpdt_calculator.getTstatTraining();
     t_testing = mdrpdt_calculator.getTstatTesting();
     mor_vals = mdrpdt_calculator.getMatchedOddsRatio();
     
     t_iter = t_stats.begin();
     
     if(t_testing.size()>0){
      for(test_iter = t_testing.begin(), m_iter=mor_vals.begin(); test_iter != t_testing.end(); t_iter++, test_iter++, m_iter++){
        currResult->analysisScores.push_back(*t_iter);
        currResult->analysisScores.push_back(*test_iter);
        currResult->analysisScores.push_back(*m_iter);
      }
     }
     else
       currResult->analysisScores.push_back(*t_iter);
       
     if(n_perms > 0)
       currResult->analysisScores.push_back(mdrpdt_calculator.getPvalue(mdrpdt_calculator.getAvgMatchedOddsRatio()));
     
     if(currResult->analysisScores.back() < 0)
      currResult->analysisScores.back() = 0;
     
    if(currResult->analysisScores.back() < threshold){
      resultList.erase(currResult++);
    }
    else 
      currResult++;
    
  } 

}
}

