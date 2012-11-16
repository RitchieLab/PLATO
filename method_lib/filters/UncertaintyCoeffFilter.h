// UncertaintyCoeffFilter.h

#ifndef __UNCERTAINTYCOEFFFILTER_H__
#define __UNCERTAINTYCOEFFFILTER_H__

#include "Filter.h"
#include <UncertaintyCoefficient.h>
#include "Permutation.h"
namespace Filters{
/// Filter calculates Uncertainty Coefficient on contingency table
class UncertaintyCoeffFilter: public Filter{

  public:
    UncertaintyCoeffFilter();
    UncertaintyCoeffFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);

    /// Set parameters for filter run
    void set_params(PARAMS params, Methods::DataSet*);

    /// Analysis for the contingency table and those filters that accept it
    void analyze(Result& res, Methods::ContingencyTable& table);

    /// returns 1 when first score is better than second one
    virtual int score_better(float first_score, float second_score){
      return first_score>second_score?1:0;}

    /// returns estimate of run time for this filter
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);

  private:
    void initialize();

    void run_permutations(Methods::DataSet& dataset);

    double get_model_fraction(int model_size);

    enum ConfigList{
      NoMatch,
      Threshold,
      UseRawScore,
      Permutations,
      TotalTest
    };

    map<string, ConfigList> ConfigMap;
    map<string, double> NullMean;
    map<string, double> NullStdDev;

    bool use_raw_score, perm_finished;
    unsigned int n_perms;
    Permutation perm;
    float threshold;
    string totalType;
    Methods::UncertaintyCoefficient uc_calculator;

};
}

#endif
