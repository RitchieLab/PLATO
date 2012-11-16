//ChiSquareFilter.h

#ifndef __MARSFILTER_H__
#define __MARSFILTER_H__

#include "Filter.h"
#include <earth.h>
#include <ContingencyTable.h>
#include "Permutation.h"
#include <time.h>

namespace Filters{

/// Filter calculates Uncertainty Coefficient on contingency table
class MarsFilter: public Filter{

  public:
    MarsFilter();
    MarsFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);

    /// Set parameters for filter run
    void set_params(PARAMS params, Methods::DataSet*);

    /// returns 1 when first score is better than second one
//    virtual int score_better(float first_score, float second_score){
//      return first_score>second_score?1:0;}

    /// returns estimate of run time for this filter
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);

  private:
    void initialize();

    void run_permutations(Methods::DataSet& dataset);

    float call_mars(int locus);

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
//    Methods::CaConChisq chi_calculator;
    Methods::earth earth_calculator;

};

}

#endif
