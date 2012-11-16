//ResultsManager.cpp
#include "ResultsManager.h"
#include <algorithm>

using namespace Methods;
namespace Filters{
///
/// Standard constructor
///
ResultsManager::ResultsManager(){
  initialize();
}

///
/// Initializes parameters for generating
/// locus combinations
/// @return
///
void ResultsManager::initialize(){
  comboIter = generator.ComboList.begin();
  outWriterSet = false;
  noMoreCombos = false;
  totalFiltersRun = 0;
  largest_size = 0;
}

///
/// Passes result list to filters after
/// creating or modifying the list of
/// results for the filter.<br>
/// The result list is passed in pieces
/// if the manager is generating any new
/// combinations of loci so that the memory
/// usage is reduced.
/// @param currFilter Filter to run
/// @param dataset DataSet to process
/// @param options ListOptions specifying how to handle loci
///        combinations
/// @return 
/// 
void ResultsManager::run_filter(Filter * currFilter, 
  DataSet & dataset, ListOptions options){

  switch(options.op){
    case PassList:
      run_analysis(currFilter, dataset, resultList);
      break;
    case RemoveFromList:
      remove_combinations(options.lociCombinations);
      run_analysis(currFilter,dataset, resultList);
      break;
    case AddToList:
      add_combinations(options.lociCombinations, currFilter, dataset); 
      break;
    default:
      throw PlatoExcept("No list operation specified in ResultsManager::run_filter");
  };
  
  totalFiltersRun++;
  
}


///
/// Runs filter on a subset of models to estimate time taken
/// to run full set of models.
/// @param currFilter Filter to test
/// @param dataset DataSet to process
/// @param options ListOptions specifying how to handle loci combinations
/// @return ProcessEstimate
///
Filter::ProcessEstimate ResultsManager::run_filter_time(Filter* currFilter,
  DataSet & dataset, ListOptions options, double num_start_models){
  
  double num_models;
 
  switch(options.op){
    case PassList:
      num_models = num_start_models;
      break;
    case RemoveFromList:
      num_models = RemoveForTimeTest(options.lociCombinations, num_start_models, dataset.num_loci());
      break;
    case AddToList:
      num_models = AddForTimeTest(options.lociCombinations, num_start_models, dataset.num_loci());
      break;
    default:
      throw PlatoExcept("No list operation specified in ResultsManager::run_filter");  
  };
  
  
  ResultSet templist;

  if(included_model_sizes.size() > 0){
    // construct temporary result list to run using largest size and some randomly 
    // sampled loci from the set
    int total_loci = dataset.num_loci();
  
    vector<int> selected;
    int randloc;
  
    largest_size = *(included_model_sizes.rbegin());
    while(selected.size() < 10){
      randloc = int(float(rand()) / RAND_MAX * total_loci);
      if((dataset.get_markers())->at(randloc)->isEnabled()){
        selected.push_back(randloc);
      }
    }
  
    
    Result newResult;
    newResult.genoCombination.resize(largest_size);
  
    // create based on largest size of model to be checked
    int startloc = 0;
    size_t currloc;
    for(startloc=0; startloc < 10; startloc++){
      currloc = startloc;
      Result newResult;
      for(int i=0; i<largest_size; i++){
        newResult.genoCombination.push_back(selected[currloc]);
        ++currloc;
        if(currloc > selected.size()-1){
          currloc = 0;
        }
      }
      templist.push_back(newResult);
    }
  
  
    ResultSet addlist;
    addlist = templist;
//  templist.push_back(templist);
    for(int i=0; i<9; i++)
      templist.insert(templist.end(), addlist.begin(), addlist.end());
  }
  
  
  Filter::ProcessEstimate estimate = currFilter->estimate_run_time(num_models, dataset, templist, largest_size);
//   estimate.seconds = estimate.seconds / 10;
  
  return estimate;
}


///
/// Returns number of models that will be included
/// @param lociCombinations size to add
/// @param num_start_models Number of models currently in list
/// @param largest_size Largest model size that has already been in list
/// @param total_loci Total loci in dataset
/// @return total number of models
///
double ResultsManager::RemoveForTimeTest(vector<unsigned int> lociCombinations, double num_start_models,
  int total_loci){

  vector<double> total_size(largest_size+1, 0.0);
  largest_size = *(included_model_sizes.rbegin());

  double total_models = num_start_models;
  
  double sum_size=0.0;
  
  // need to remove some portion of the models based on 
  // the sizes removed
  for(int i=largest_size; i >= 1; i--){
    total_size[i] = generator.calc_combinations(total_loci, i);
    sum_size += total_size[i];
  }

  for(size_t i=0; i < lociCombinations.size(); i++){
    included_model_sizes.erase(lociCombinations[i]);
  }
  
  largest_size = *(included_model_sizes.rbegin());
  // remove the indicated percentage of sizes
  for(size_t i=total_size.size()-1; i > size_t(largest_size); --i){
    total_models -= (total_size[i] / sum_size * num_start_models);
  }

  return total_models;
}


///
/// Returns number of models that will be included 
/// @param lociCombinations size to add
/// @param num_start_models Number of models currently in list
/// @param largest_size Largest model size that has already been in list
/// @param total_loci Total loci in dataset
/// @return total number of models
///
double ResultsManager::AddForTimeTest(vector<unsigned int> lociCombinations, double num_start_models,
  int total_loci){

  for(size_t i=0; i < lociCombinations.size(); i++){
    included_model_sizes.insert(lociCombinations[i]);
  }  
  
  double total_models = num_start_models;
  int loci_in_set;
  // when empty simply add model totals  
  if(num_start_models == 0){
    loci_in_set = total_loci;
  }
  else{
    // when already existing exist use the models that have passed through
    // estimate number of unique SNPS
    if(largest_size == 1){
      loci_in_set = int(num_start_models);
    }
    else{ // figure out how many possible SNPS
      
      if(num_start_models * largest_size > total_loci){
        loci_in_set = total_loci;
      }
      else{
        loci_in_set = int(total_loci / (num_start_models * largest_size) * total_loci);
      }
    }
  }
  
  for(size_t i=0; i < lociCombinations.size(); i++){
    total_models += generator.calc_combinations(loci_in_set, lociCombinations[i]);
  }
  
  return total_models;
}



/// 
/// Used to list loci combinations in current list
///
void ResultsManager::list_combinations(ResultSet & resList){
  for(ResultIter currResult = resList.begin(); currResult !=resList.end();){
    for(unsigned int x=0; x<(*currResult).genoCombination.size(); x++){
      cout << (*currResult).genoCombination[x] << " ";
    }
    cout << endl;
    ++currResult;
  }
  cout << endl;
}

///
/// Adds combinations of loci to the list<br> 
/// Combinations are based on the presence of loci in the
/// current set.  If no loci are present all the loci 
/// are used.
/// @param lociCombinations vector Size of combinations to add
/// @param currFilter Filter for analysis
/// @param dataset DataSet for this analysis
/// @return
///
void ResultsManager::add_combinations(vector<unsigned int> & lociCombinations,
  Filter * currFilter, DataSet & dataset){
  vector<unsigned int> lociInSet;
  unsigned int numLoci;
  bool ran_previously = false;
  // set the loci for use in adding combinations
  // when list is empty and this is first filter use all the loci 
  if(resultList.empty() && totalFiltersRun==0){
     numLoci = dataset.num_loci(); 
     for(unsigned int currLoc=0; currLoc < numLoci; currLoc++){
       lociInSet.push_back(currLoc);
     }
  }
  else{
    set_included_loci(lociInSet); 
  }

  // use lociInSet to generate combinations for addition to 
  // the list 
  generator.SetLoci(lociInSet.size());  

  generator.SetComboInterval(DEFAULTCOMBOINTERVAL);   

  unsigned int numCombos = lociCombinations.size();

  // transfer results list to a temporary list
  ResultSet hold_results;
  hold_results.splice(hold_results.begin(), resultList); 
  ResultSet filter_arguments;
  ResultIter holdIter = hold_results.begin();
  ResultIter transfer_start, transfer_end;

  Result newResult;
  if(totalFiltersRun)
    newResult.analysisScores.assign(totalFiltersRun, -1e30);
  for(unsigned int currCombSize=0; currCombSize<numCombos; currCombSize++){
    // check that you have enough loci in list to do this combination
    // else skip to next
    if(lociInSet.size() < lociCombinations[currCombSize])
      continue;
    
    noMoreCombos = false;
    generator.ComboEnds(lociCombinations[currCombSize], lociCombinations[currCombSize]);
    get_next_combination(newResult.genoCombination, lociInSet); // always have at least one combination
    transfer_start = hold_results.begin();
    transfer_end = hold_results.end();
    for(holdIter = transfer_start; holdIter != hold_results.end(); holdIter++){
      if(holdIter->genoCombination.size() >= lociCombinations[currCombSize]){
        transfer_end = holdIter;
        break;
      }
    }
    
    // transfer any results that were found
    if(transfer_start != hold_results.end()){
      filter_arguments.splice(filter_arguments.begin(), hold_results, 
        transfer_start, transfer_end);  
    }
    // fill argument list with as many arguments as needed up to the 
    // limit or all combinations of this size are used
    do{
      if(holdIter != hold_results.end()){
        if(holdIter->genoCombination != newResult.genoCombination){
          filter_arguments.push_back(newResult);
        }
        else{
          filter_arguments.splice(filter_arguments.end(), hold_results, holdIter);
        }
      }
      else{
        filter_arguments.push_back(newResult);    
      } 
      if(filter_arguments.size() >= DEFAULTCOMBOINTERVAL){
        run_analysis(currFilter, dataset, filter_arguments);
        ran_previously=true;
        // splice the filter_arguments onto the resultList
        resultList.splice(resultList.end(), filter_arguments); 
      }
    }while(get_next_combination(newResult.genoCombination, lociInSet)); 
  } 

  // if any left on hold_results transfer to filter_arguments and
  // run if necessary  
  filter_arguments.splice(filter_arguments.end(), hold_results); 
  

  if(!filter_arguments.empty()){
    run_analysis(currFilter, dataset, filter_arguments);
    resultList.splice(resultList.end(), filter_arguments);
  }
  else if(!ran_previously){
  }
  
}


///
/// Returns next combination from combination generator.<br>
/// Continues to generate combinations as requested <br>
/// Returns false when all combinations have been returned.
/// @param newCombination LocusVec that is set to the new combination
/// @param lociInSet LocusVec listing all loci in the current set
///        so that loci can be transformed back into original IDs
/// @return returns true when a new combination is set and false
/// when all combinations have been generated
///
bool ResultsManager::get_next_combination(LocusVec & newCombination,
  LocusVec & lociInSet){
  if(comboIter != generator.ComboList.end()){
    unsigned int comboSize = (*comboIter).size();
    for(unsigned int currLoc=0; currLoc < comboSize; currLoc++)
      (*comboIter)[currLoc] = lociInSet[(*comboIter)[currLoc]];
    newCombination = *comboIter;
    comboIter++;
    return true; 
  }
  
  if(!noMoreCombos){
//if(generator.param_AlreadyStarted()) logger.log_parameters(generator, "plato_combo.log");
    noMoreCombos = generator.GenerateCombinations();
    comboIter = generator.ComboList.begin();
  }
  else{
    return !noMoreCombos;
  }

  if(!generator.ComboList.empty()){
    unsigned int comboSize = (*comboIter).size();
    for(unsigned int currLoc=0; currLoc < comboSize; currLoc++)
      (*comboIter)[currLoc] = lociInSet[(*comboIter)[currLoc]];

    newCombination = *comboIter;
    comboIter++;
    return true;
  }
  else{
    return false;
  }
}


///
/// Iterates through results and stores loci 
/// @param lociInSet vector that will contain indexes into locus set
/// @return
///
void ResultsManager::set_included_loci(LocusVec & lociInSet){
  map<int,int> includedLoci;
  map<int,int>::iterator mapIter;
  ResultIter resIter, endIter = resultList.end();
  vector<unsigned int>::iterator locIter;
  for(resIter = resultList.begin(); resIter != endIter;resIter++){
    for(locIter = resIter->genoCombination.begin();
      locIter != resIter->genoCombination.end(); locIter++){
      includedLoci[*locIter]=1;
    }
  } 
  
  for(mapIter = includedLoci.begin(); mapIter != includedLoci.end();
    mapIter++){
    lociInSet.push_back(mapIter->first);  
  }

  // sort the loci 
  sort(lociInSet.begin(), lociInSet.end(), less<int>());
  
}


///
/// Runs actual analysis
/// @param currFilter Filter to use 
/// @param dataset DataSet for this analysis
/// @param results ResultSet list for storing analysis results
/// @return
///
void ResultsManager::run_analysis(Filter * currFilter, DataSet & dataset,
  ResultSet & results){
  currFilter->analyze(results, dataset); 
}


///
/// Removes all the results in the list matching the
/// sizes specified
/// @param lociCombinations vector containing list of combination
/// sizes to remove from the list before analyzing
/// @return
///
void ResultsManager::remove_combinations(vector<unsigned int> & lociCombinations){

  // list is ordered from smallest size to largest so need
  // to find the start and finish of each and then remove
  // them from list 
  unsigned int numSizes = lociCombinations.size();
  bool startMarked = false;
  ResultIter resIter, startErase, endErase, endList;
 
  for(unsigned int startSize=0; startSize < numSizes; startSize++){
    endList = resultList.end();
    startMarked = false;
    for(resIter = resultList.begin(); resIter != endList; resIter++){
      if(resIter->genoCombination.size() == lociCombinations[startSize]){
        if(startMarked) 
          endErase = resIter;
        else{
          startErase = endErase = resIter;
          startMarked = true;
        }
      } 
      else if(resIter->genoCombination.size() > lociCombinations[startSize]){
        if(startMarked){
          resultList.erase(startErase, ++endErase);
          startMarked = false;
        }
        break;
      }
    } 
    if(startMarked)
      resultList.erase(startErase, ++endErase);
  } 
}
}

