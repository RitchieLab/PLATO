//MaxScoreFilter.cpp

#include "MaxScoreFilter.h"
#include <iostream>
#include <Helper.h>

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
MaxScoreFilter::MaxScoreFilter():Filter("MaxScore Filter"){
  initialize();
}

///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
MaxScoreFilter::MaxScoreFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void MaxScoreFilter::initialize(){
  n_perms = 0;
  perm_completed = false;
  threshold=1.01;
  ConfigMap["PERMUTATIONS"] = Permutations;
  ConfigMap["THRESHOLD"] = Threshold;
  name_adjusted = false;
}


///
/// Calculates the p values for each subfilter and retains the 
/// lowest p value which is the maximum score.  If permutation testing
/// is also included will add a score based on the permutation score
/// for that p value.
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void MaxScoreFilter::analyze(ResultSet & resultList, DataSet & dataset){

   unsigned int curr_sub;
   unsigned int num_subs = subfilters.size();

   set_text_outname();

  // first establish permutation distributions for all the subfilters
  // on first set -- don't do it on subsequent sets
  if(!perm_completed){
    run_perms(dataset);
    set_text_outname();
  }
  
  for(curr_sub=0; curr_sub < num_subs; curr_sub++){
      subfilters[curr_sub]->analyze(resultList, dataset);
  }
  
  
  multiset<float>::iterator distIter;
  vector<float>::iterator scoreIter;
  
  // second for each result in list -- analyze and set permutation derived
  // p value 
  float min;
  for(ResultIter resIter = resultList.begin(); resIter !=resultList.end();){ 
    min = 1.00;
    for(scoreIter = resIter->analysisScores.begin(); scoreIter != resIter->analysisScores.end();
      scoreIter++){
      if(*scoreIter < min){
        min = *scoreIter;
      }
    }

    resIter->analysisScores.clear();
    
    if(min < threshold){
      
      resIter->analysisScores.push_back(min);
    
      // determine p value of distribution
      if(!distributions.empty()){
        unsigned int better=0;
        for(distIter = distributions[resIter->genoCombination[0]].begin(); 
          distIter!=distributions[resIter->genoCombination[0]].end(); distIter++){
          if(*distIter < min){
            better++;
          }
        }
        resIter->analysisScores.push_back(float(better)/distributions[resIter->genoCombination[0]].size());
      }
      resIter++;
    }
    else{
      resultList.erase(resIter++);
    }
  }
 
}


///
/// Run permutations and create distributions
///
void MaxScoreFilter::run_perms(DataSet & data){

  // store original status
  vector<bool> original_status;
  vector<unsigned int> all_inds;
  int num_inds = data.num_inds();
  for(int ind=0; ind < num_inds; ind++){
    original_status.push_back(data.get_sample(ind)->getAffected());
    all_inds.push_back(ind);
  } 
  
  unsigned int num_affected = data.num_affected();
 
  multiset<float> sampleset;
  // create new multiset vector for p value distributions
//   distributions.assign(multiset<float>, data.num_loci());
  distributions.assign(data.num_loci(),sampleset);
 
  ResultSet permList;
  Result temp;
  temp.genoCombination.resize(1);
  
  unsigned int num_subs = subfilters.size(), curr_sub;
  
  for(unsigned int curr_perm=0; curr_perm < n_perms; curr_perm++){
    perm.permute_status(data, all_inds, num_affected);
// perm.output_set(data, curr_perm);    
    permList.clear();
    
//     reorderAlleles(samples, markers);
    reorderAlleles(data.get_samples(), data.get_markers());
    
    
    // construct list for analysis
    for(unsigned int currloc=0; currloc < data.num_loci(); currloc++){
      temp.genoCombination[0] = currloc;
      permList.push_back(temp);
    }
    
    // need to run and save all the results
    for(curr_sub=0; curr_sub < num_subs; curr_sub++){
      subfilters[curr_sub]->analyze(permList, data);      
    }
    
    // find best for each and store in appropriate 
    vector<float>::iterator scoreIter;
    unsigned int currloc=0;
    
    for(ResultIter resIter = permList.begin(); resIter != permList.end(); resIter++){
      float min = 1.0;
      for(scoreIter = resIter->analysisScores.begin(); scoreIter != resIter->analysisScores.end();
        scoreIter++){
        if(*scoreIter < min){
          min = *scoreIter;
        }
      }
      distributions[currloc++].insert(min);
    }
  }
  
  
  // return status to original values
  perm.set_status(original_status, data);
  reorderAlleles(data.get_samples(), data.get_markers());
}



///
/// Sets name of this filter to reflect sub filters
/// used -- Will appear as column headers in text file
/// @return
///
void MaxScoreFilter::set_text_outname(){
//   if(subfilters.size() == 0)
//     return;
//   name = subfilters[0]->getName();
//   unsigned int num_subs = subfilters.size();
//   for(unsigned int curr_sub=1; curr_sub < num_subs; curr_sub++){
//       name += "\t" + subfilters[curr_sub]->getName();
//   }
  
  if(!name_adjusted && n_perms >0){
    name += "\tPermuted p-value";
    name_adjusted = true;
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
Filter::ProcessEstimate MaxScoreFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  
  Filter::ProcessEstimate time_estimate;

  clock_t start, end;
  

  size_t num_run = resultList.size();
 
  unsigned int original_perms = n_perms;
  n_perms = 0;
  start = clock();
  analyze(resultList, dataset);
  end = clock();
  n_perms = original_perms;
  perm_completed = false;

  time_estimate.seconds = double(end - start) / (CLOCKS_PER_SEC * num_run) ;
  time_estimate.seconds *= num_models * (n_perms+1);
  
    time_estimate.num_models = threshold * num_models;
     
  return time_estimate;
}


///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void MaxScoreFilter::set_params(PARAMS params, DataSet* set){
  std::map<string, string>::iterator configIter;
  string report_loci_str;
  stringstream ss;
//   unsigned int locus;
  
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

