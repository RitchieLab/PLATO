// ResultsManager.h

#ifndef __RESULTSMANAGER_H__
#define __RESULTSMANAGER_H__

#include "Filter.h"
#include <ComboGenerator.h>
#include <filters/ResultSet.h>
#include <filters/OutputWriter.h>
#include <filters/PlatoDefinitions.h>
#include <filters/LogInfo.h>
#include <set>

namespace Filters{

#define DEFAULTCOMBOINTERVAL 20000

///
/// Manages the results for PLATO<br>
/// Can add or remove combinations from the list<br>
/// Passes previous results to filters for processing
/// in segments so no need to generate all combinations
/// of loci at one time
/// 

/// Maintains list of results for plato for passing to filters
class ResultsManager{

  public:
  
    /// Standard constructor
    ResultsManager();
    
    /// destructor
    virtual ~ResultsManager(){}
    
    /// passes results to filters for processing
    virtual void run_filter(Filter * currFilter, Methods::DataSet & dataset,
      ListOptions currentOp);

    /// Lists combinations to standard out
    void list_combinations(ResultSet & resList);

    ResultSet resultList;
    
    void set_proc_num(int num){proc_num=num;}
    int get_proc_num(){return proc_num;}

    Filter::ProcessEstimate run_filter_time(Filter* currFilter,
      Methods::DataSet & dataset, ListOptions options, double num_start_models);

  protected:
    int proc_num, largest_size;
    set<int> included_model_sizes;

   void run_analysis(Filter * currFilter, Methods::DataSet & dataset,
	    ResultSet & results);
    void remove_combinations(vector<unsigned int> & lociCombinations);
    void add_combinations(vector<unsigned int> & lociCombinations, Filter * currFilter, 
      Methods::DataSet & dataset);
    void set_included_loci(LocusVec & lociInSet);
    void initialize();
    bool get_next_combination(LocusVec & newCombination, LocusVec & locInSet);


    double RemoveForTimeTest(vector<unsigned int> lociCombinations, double num_start_models,
      int total_loci);

    double AddForTimeTest(vector<unsigned int> lociCombinations, double num_start_models,
      int total_loci);

    Methods::ComboGenerator generator;
    vector<LocusVec>::iterator comboIter;
    OutputWriter * outWriter;
    LogInfo logger;
    unsigned int totalFiltersRun;
    bool outWriterSet, noMoreCombos;
};
}

#endif
