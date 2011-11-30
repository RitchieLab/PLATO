// Permutation.cpp

#include "Permutation.h"

#include <iostream>

using namespace Methods;
namespace Filters{
// #include <time.h>

///
/// Constructor
///
Permutation::Permutation(){
  initialize();
}

///
///  Sets the number of permutations
///  @param n_perms Number of permutations to run
///  @return 
///
void Permutation::set_num_permutations(int n_perms){
  num_perms = n_perms;
}

///    
/// Adds filter to list for running permutation test
/// @param permuted_filter Filter to add to establish permuted 
/// distribution
/// @return 
///
void Permutation::add_filter(Filter * permuted_filter){
  filters.push_back(permuted_filter);
  sort_priorities.push_back(permuted_filter->get_sort_priority());
  p_value_distribution.push_back(vector<float>(0,0.0));
}
   

///   
/// Returns p value for indicated value
/// Assumes only one filter has been set
/// @param value Value from filter 
/// @param filter_number unsigned int containing filter number
/// to check against -- 0 is default
/// @return pv value
///
float Permutation::get_p_value(float value, unsigned int filter_number){

  float p_value = 0.0;
  int index=0;
  
  for(;index < num_perms; index++){
    if(value < p_value_distribution[filter_number][index]){
      break;
    }
  }
  
  switch(sort_priorities[filter_number]){
    case Filter::highest:
      p_value = (num_perms-index)/float(num_perms);
      break;
    case Filter::lowest:
      p_value = index/float(num_perms);
      break;
  }
  
  return p_value;
}

/// 
/// Runs permutations and creates distribution for contingency 
/// @param data DataSet containing set to permute
/// @return 
///
void Permutation::run_permutations_contingency(DataSet & data){
  
  // store original status vectors
//   vector<unsigned int> original_affected = data.get_affected_vector();
//   vector<unsigned int> original_unaffected = data.get_unaffected_vector();
//   unsigned int num_affected = original_affected.size();

  // store original status
  vector<bool> original_status;
  vector<unsigned int> all_inds;
  int num_inds = data.num_inds();
  for(int ind=0; ind < num_inds; ind++){
    original_status.push_back(data.get_sample(ind)->getAffected());
    all_inds.push_back(ind);
  } 

  unsigned int num_affected = data.num_affected();
 
  for(int curr_perm=0; curr_perm < num_perms; curr_perm++){
    permute_status(data, all_inds, num_affected);
    get_best_scores(data);
if((curr_perm+1) % 10 == 0)
cout << "Completed permutation " << curr_perm +1 << endl;
  }

  // return status to original values
  set_status(original_status, data);

  // sort the distributions in ascending order
  for(unsigned int distribution=0; distribution<p_value_distribution.size();
    distribution++){
    sort(p_value_distribution[distribution].begin(), p_value_distribution[distribution].end());
  }

}


///
/// Sets the best scores for each of the filters
/// set within this class
/// @param data permuted DataSet
/// @return
///
void Permutation::get_best_scores(DataSet & data){
  
  vector<Filter *>::iterator filterIter;
  
  vector<float>  best_scores;
  for(unsigned int i=0; i<filters.size(); i++){
    best_scores.push_back( filters[i]->get_worst_score());
  }


  // currently only works on single loci
  int num_loci = data.num_loci();
  
  ContingencyTable table;
  unsigned int currFilter;
  unsigned int num_filters = filters.size();
  
  for(int curr_loc=0; curr_loc < num_loci; curr_loc++){
    table.get_counts(curr_loc, &data);
    
    Result newResult;
    newResult.genoCombination.push_back(curr_loc);

    for(currFilter=0; currFilter < num_filters; currFilter++){   
      filters[currFilter]->analyze(newResult, table);
      if(filters[currFilter]->score_better(newResult.analysisScores[currFilter], 
        best_scores[currFilter])){
          best_scores[currFilter] = newResult.analysisScores[currFilter];
        }
    }

  }
    
  for(unsigned int i=0; i<p_value_distribution.size(); i++){
    p_value_distribution[i].push_back(best_scores[i]);
  }
  
  
  
}


void Permutation::output_set(DataSet& data, int perm_num){
  cout << "perm start "<< perm_num << endl;
  
  Sample * ind;
  
  for(int i =0; i<data.num_inds(); i++){
    ind = data[i];
    cout << int(ind->getAffected());
    for(unsigned int currloc=0; currloc < data.num_loci(); currloc++){
      cout << " " << ind->get_genotype(currloc);
    }
    cout << endl;
  }
  cout << endl;
}


///    
/// Runs permutations and creates distribution
/// @param data DataSet containing set to permute
/// @return
///
void Permutation::run_permutations(DataSet& data){
  
  // store original status vectors
//   vector<unsigned int> original_affected = data.get_affected_vector();
//   vector<unsigned int> original_unaffected = data.get_unaffected_vector();
//   unsigned int num_affected = original_affected.size();
  
  
  // store info on affected and unaffected
//   vector<unsigned int> original_affected;
//   vector<unsigned int> original_unaffected;
//   vector<unsigned int> all_inds;
//   int num_inds = data.num_inds();
//   for(int ind=0; ind < num_inds; ind++){
//     if(data.get_sample(ind)->getAffected())
//       original_affected.push_back(ind);
//     else
//       original_unaffected.push_back(ind);
//     all_inds.push_back(ind);
//   }
  
  // store original status
  vector<bool> original_status;
  vector<unsigned int> all_inds;
  int num_inds = data.num_inds();
  for(int ind=0; ind < num_inds; ind++){
    original_status.push_back(data.get_sample(ind)->getAffected());
    all_inds.push_back(ind);
  }  


  unsigned int num_affected = data.num_affected();
  
  ResultSet perm_results;
  Result curr_result;
  curr_result.genoCombination.push_back(0);
  
  // currently only works on single loci
  int num_loci = data.num_loci();
  int total;
  
  float best_score;
  
  // run permutations to create distribution
  for(int curr_perm=0; curr_perm < num_perms; curr_perm++){
    // permute dataset
    permute_status(data, all_inds, num_affected);

    best_score = filters[0]->get_worst_score();
    // currently only works on single locus
    total=0;
    
    ComboGenerator generator;
    generator.ComboEnds(min_model_size, max_model_size);
    generator.SetLoci(num_loci);  
    generator.SetComboInterval(30000);   
    
    bool done = generator.GenerateCombinations();
    std::vector < std::vector<unsigned int> >::iterator combo_iter;
    while(generator.ComboList.size() > 0){
      for(combo_iter = generator.ComboList.begin(); combo_iter != generator.ComboList.end();
        combo_iter++){
        curr_result.genoCombination = *combo_iter;
        perm_results.push_back(curr_result);
      }
    
      generator.ComboList.clear();
      
      // run results
      filters[0]->analyze(perm_results, data);    
      // check results and keep best
      for(ResultIter currResult = perm_results.begin(); currResult != perm_results.end(); currResult++){ 
        if(filters[0]->score_better(currResult->analysisScores[0], 
          best_score)){
            best_score = currResult->analysisScores[0];
        }
      }
      
      perm_results.clear();
      if(!done)
        generator.GenerateCombinations(); 
    } 
    p_value_distribution[0].push_back(best_score);
    
  }
  
  sort(p_value_distribution[0].begin(), p_value_distribution[0].end());

  // return status to original values
  set_status(original_status, data);

}



///    
/// Runs indicated number of permutations and creates distribution
/// @param data DataSet containing set to permute
/// @param n_perms Number of permutations to run
/// @return
///
void Permutation::run_permutations(DataSet& data, unsigned int n_perms){
  num_perms = n_perms;
  run_permutations(data);
}

///
/// Initializes variables
///
void Permutation::initialize(){
  num_perms = 1000;

  min_model_size = 1; 
  max_model_size = 1;

  elapsed_time=0.0;

}
    

///
/// Sets status according to vector passed
/// @param status
/// @param data
///
void Permutation::set_status(vector<bool>& status, DataSet& data){
  int num_inds = data.num_inds();
  for(int ind=0; ind < num_inds; ind++){
    data.get_sample(ind)->setAffected(status[ind]);
    data.get_sample(ind)->setPheno(status[ind]+1);
  }   
}


///
/// Permutes status 
/// @param data DataSet to permute
/// @param num_affected unsigned int Number of affected to set
/// @return 
///
void Permutation::permute_status(DataSet& data, vector<unsigned int>& all_inds,
  unsigned int num_affected){
  
  vector<unsigned int> affected, unaffected;
  unsigned int last_ind = data.num_inds()-1;

  vector<bool> status_vec(data.num_inds(), false);

  unsigned int selected_index, temp, last_active=0, curr_ind;
  
  for(curr_ind=0; curr_ind < num_affected; curr_ind++){
    last_active = last_ind-curr_ind;
    selected_index = rand() % last_active;
    temp = all_inds[selected_index];
    all_inds[selected_index] = all_inds[last_active];
    all_inds[last_active] = temp;
    status_vec[temp] = true;
  }
  
  set_status(status_vec, data);
  

}
}
 
