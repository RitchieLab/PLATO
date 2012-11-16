// LogRegFilter.h

#ifndef __LOGREGRESSFILTER_H__
#define __LOGREGRESSFILTER_H__

#include "Filter.h"
#include <FlatIndex.h>
#include <LogisticRegression.h>

///
/// Filter conducts logistic regression.
/// Can run on multiple locus combinations
///

namespace Filters{

/// Runs logistic regression
class LogRegFilter: public Filter{
  
  public:
    LogRegFilter();
    LogRegFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);

    /// returns estimate of run time for this filter
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
    
  private:
    void initialize();

    enum ConfigList{ 
      NoMatch,
      Threshold,
//       CalcType,
      ModelType,
//       ComboLimit,
//       ComboMin,
      MaximumIterations,
      InteractionIncluded,
      PType,
      IncludeCoeff,
      IncludeIntercept,
      MaxModelSize
    };
    
    // specifies which p value to use
    enum PValueTypes{
      Overall,
      Coefficient,
      None
    };
    
    float pThreshold;
    PValueTypes pValType;
    map<string, ConfigList> ConfigMap;
    map<string, int> TotalTypeDF;
    string modType;
    map<string, PValueTypes> PValueTypeMap;
    Methods::LogisticRegression lr_calculator;
    bool coeff_output, coeff_intercept_out, interactIncluded;
    int max_model_size, max_coeff;
    
};
}
#endif
