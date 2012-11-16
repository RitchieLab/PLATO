 // Filter.h


#ifndef __FILTER_H__
#define __FILTER_H__

#include <map>
#include <list>
#include <DataSet.h>
#include <filters/Stringmanip.h>
#include "FilterExcept.h"
#include "ResultSet.h"
#include "FilterFactory.h"
#include "FilterDefinitions.h"
#include <ContingencyTable.h>


///
/// Abstract base class for different filters used in PLATO project
/// Current filter names: <br>
/// CHISQUARE <br>
/// SERF <br>
/// MDR <br>
/// LOGISTICREGRESS
///

namespace Filters{

/// Abstract base class for all analysis
class Filter{

  public:  
    Filter(string filterName);
    Filter(const Filter & origFilter);
    virtual ~Filter();
    virtual Filter & operator=(const Filter & origFilter);
    
    ///
    /// Must be defined in each child class. 
    /// Analyzes data using filter and stores results in resultList
    /// Loci or combinations of loci that do not meet criteria for inclusion
    /// can be deleted from the resultList
    /// @param resultList ResultSet for storing results of analysis
    /// @param dataset DataSet containing data for analysis
    /// @return none
    ///
    virtual void analyze(ResultSet & resultList, Methods::DataSet & dataset)=0;
    
    ///
    /// Must be defined in each child class. 
    /// Sets configuration parameters for the filter
    /// Should throw an exception when parameter is invalid
    /// @param params Map with configuration keyword and value as string
    /// @param set DataSet
    /// @return none
    /// @throws FilterExcept when any parameters passed are invalid
    ///
    virtual void set_params(PARAMS params, Methods::DataSet* set)=0;
    
    
    /// Returns name of this filter
    inline string getName() const {return name;}
    
    /// adds subfilter to this filter
    virtual void add_subfilter(Filter * subfilter);
    
    /// sets the analysisID to the filter for use in logging 
    /// results for filters that implement it
    virtual void set_filterID(string anID){filterID=anID;}

    /// Analysis for the contingency table and those filters that accept it
    virtual void analyze(Result& res, Methods::ContingencyTable& table){}
    
    /// returns starting score for use in permutation testing comparisons
    virtual float get_worst_score(){return worst_score;}
    
    /// returns 1 when first score is better than second one (less than in this case)
    virtual int score_better(float first_score, float second_score){
      return first_score<second_score?1:0;}
    
    /// Enumeration indicating best score is lowest or highest
    enum BestScore{
      lowest,
      highest
    };
    
    /// Returns whether lower or higher results are better
    /// for raw scores generated in filter
    BestScore get_sort_priority(){return sort_priority;}
    
    struct ProcessEstimate{
      double seconds;
      double num_models;
    };
    
    //Returns number ProcessEstimate with number of models expected and time to evaluate
    virtual ProcessEstimate estimate_run_time(double num_models,Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1){
        ProcessEstimate est; 
        est.seconds=0; 
        est.num_models=num_models;
        return est;
     }
    
  protected:
  
    double normdist(double xx, double mean, double standard_dev);
//     double timespecdiff(timespec* x, timespec* y);
  
    string name, filterID;
    vector<Filter *> subfilters;
    BestScore sort_priority;
    float worst_score;
  
  private:
    void copy(const Filter & origFilter);
    
};

}

#endif
