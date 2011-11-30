//ContingencyFilter.cpp

#include "ContingencyFilter.h"

using namespace Methods;


namespace Filters{
///
/// constructor -- initialize variables
///
ContingencyFilter::ContingencyFilter():Filter("Contingency Filter"){
  initialize();
}

///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
ContingencyFilter::ContingencyFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void ContingencyFilter::initialize(){
  n_perms = 1000;
  perm_completed = false;
  threshold=1.01;
  ConfigMap["PERMUTATIONS"] = Permutations;
  ConfigMap["THRESHOLD"] = Threshold;
  ConfigMap["ALWAYSINCLUDE"] = ReportLoci;
}


///
/// Calculates contingency tables for all
/// subfilters that use the contingency table class.
/// Optimized for quicker speed.  Currently can
/// only be used with the text output.
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void ContingencyFilter::analyze(ResultSet & resultList, DataSet & dataset){

   unsigned int curr_sub;
   unsigned int num_subs = subfilters.size();

set_text_outname();

   // first establish permutation distributions for all the subfilters
   // on first set -- don't do it on subsequent sets
  if(!perm_completed){
    for(curr_sub=0; curr_sub < num_subs; curr_sub++){
      perm.add_filter(subfilters[curr_sub]);
    }
    perm.set_num_permutations(n_perms);
    perm.run_permutations_contingency(dataset);
    perm_completed = true;
    set_text_outname();
  }

  ContingencyTable table;
  
  vector<Marker*>* markers = dataset.get_markers();
  int locus;

  // second for each result in list -- analyze and set permutation derived
  // p value 
  for(ResultIter currResult = resultList.begin(); currResult !=resultList.end(); ++currResult){ 
    locus = (*markers)[currResult->genoCombination[0]]->getLoc();
    table.get_counts(locus, &dataset);
//     table.get_counts(currResult->genoCombination[0], &dataset);
    
    for(curr_sub=0; curr_sub < num_subs; curr_sub++){
      subfilters[curr_sub]->analyze(*currResult, table);
      currResult->analysisScores.back() = perm.get_p_value(currResult->analysisScores.back(), curr_sub);
    }
  }
  
  unsigned int start_index = resultList.begin()->analysisScores.size() -  num_subs;
  unsigned int end_index = resultList.begin()->analysisScores.size();
  bool keep_result;
  
  // Finally compare the p values to the specified threshold
  // and output only those loci who have at least one result better than
  // the threshold -- set others to -1 
  // Also keep those who are found in the map for reporting certain loci
  for(ResultIter currResult = resultList.begin(); currResult != resultList.end();){
    keep_result = false;
    
    if(report_loci.find(currResult->genoCombination[0]) != report_loci.end())
      keep_result = true;
    else
      for(curr_sub = start_index; curr_sub < end_index; ++curr_sub){
        if(currResult->analysisScores[curr_sub] <= threshold)
          keep_result = true;
        else
          currResult->analysisScores[curr_sub] = -1;
      }
      
    if(!keep_result)
      resultList.erase(currResult++);
    else
      ++currResult;
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
Filter::ProcessEstimate ContingencyFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  

  Filter::ProcessEstimate sub_estimate, time_estimate;
  sub_estimate.num_models = num_models;
  double time_used = 0.0;
  
  // ResultSet testList = resultList;
  
  // need to call estimate on each sub filter
  for(size_t curr_sub=0; curr_sub < subfilters.size(); curr_sub++){
    ResultSet subList = resultList;
    sub_estimate = subfilters[curr_sub]->estimate_run_time(sub_estimate.num_models, 
      dataset, subList, 1);
    time_used += sub_estimate.seconds;
  }

  time_estimate.seconds = time_used;
  time_estimate.num_models = sub_estimate.num_models;
  
  return time_estimate;
  
}



///
/// Sets name of this filter to reflect sub filters
/// used -- Will appear as column headers in text file
/// @return
///
void ContingencyFilter::set_text_outname(){
  if(subfilters.size() == 0)
    return;
  name = subfilters[0]->getName();
  unsigned int num_subs = subfilters.size();
  for(unsigned int curr_sub=1; curr_sub < num_subs; curr_sub++){
      name += "\t" + subfilters[curr_sub]->getName();
  }  
}


///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void ContingencyFilter::set_params(PARAMS params, DataSet* set){
  std::map<string, string>::iterator configIter;
  string report_loci_str;
  stringstream ss;
  unsigned int locus;
  
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
      case ReportLoci:
        report_loci_str = configIter->second;
        ss.str(report_loci_str);
        while(!ss.eof()){
          ss >> locus;
          locus--;
          report_loci[locus]=true;
        }
        break;
      case Permutations:
        n_perms = Stringmanip::stouint(configIter->second);
        if(n_perms < 1 || n_perms > 100000)
          throw FilterExcept(configIter->first + " must be between 1 and 100000 for " +
            name + " filter\n\n");
        break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  }  
}
}

