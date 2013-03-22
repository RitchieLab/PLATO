//ArmitageFilter.h

#ifndef __ARMITAGEFILTER_H__
#define __ARMITAGEFILTER_H__

#include "Filter.h"
#include "../ChiSquareArmitage.h"
#include <ContingencyTable.h>
#include "Permutation.h"
#include <time.h>

namespace Filters{

/// Filter calculates Uncertainty Coefficient on contingency table
class ArmitageFilter: public Filter{

  public:
    ArmitageFilter();
    ArmitageFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    
    /// Set parameters for filter run
    void set_params(PARAMS params, Methods::DataSet*);
  
    /// Analysis for the contingency table and those filters that accept it
    void analyze(Result& res, Methods::ContingencyTable& table);
  
    /// returns 1 when first score is better than second one
    virtual int score_better(float first_score, float second_score){
      return first_score>second_score?1:0;}

    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
  
  private:
    void initialize();
    
    void run_permutations(Methods::DataSet& dataset);
    
    float call_armitage(Methods::ContingencyTable& table);
    
    enum ConfigList{ 
      NoMatch,
      Threshold,
      UseRawScore,
      Permutations,
      TotalTest
    };
    
    map<string, ConfigList> ConfigMap;
    map<string, Methods::ContingencyTable::TotalType> TotalTypeMap;
    
    bool use_raw_score, perm_finished;
    unsigned int n_perms;
    Permutation perm;
    float threshold;
    Methods::ContingencyTable::TotalType total_type;
    
    float arm_score;
    Methods::ChiSquareArmitage arm_calculator;

};

}
#endif
