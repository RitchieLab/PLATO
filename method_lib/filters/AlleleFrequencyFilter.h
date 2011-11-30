//AlleleFrequencyFilter.h

#ifndef __ALLELEFREQFILTER_H__
#define __ALLELEFREQFILTER_H__

#include "Filter.h"
#include <AlleleFrequency.h>
#include <time.h>

///
/// Filter conducts allele frequency calculation.
///

namespace Filters{

/// Runs logistic regression
class AlleleFrequencyFilter: public Filter{
  
  public:
    AlleleFrequencyFilter();
    AlleleFrequencyFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);

    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);

  private:
    void initialize();

    enum ConfigList{ 
      NoMatch,
      Threshold,
      IncludeType
    };

    float freqThreshold;
    Methods::AlleleFrequency freq_calculator;
    map<string, ConfigList> ConfigMap;
    map<string, string> FreqMap;
    string frequency_option;

};

}

#endif
