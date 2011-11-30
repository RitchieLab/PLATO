// MDRFilter.h

#ifndef __MDRFILTER_H__
#define __MDRFILTER_H__

#include "Filter.h"
#include <FlatIndex.h>
#include <MDR.h>
#include <time.h>

///
/// Filter conducts logistic regression.
/// Can run on multiple locus combinations
///
namespace Filters{
/// Runs logistic regression
class MDRFilter: public Filter{
  
  public:
    MDRFilter();
    MDRFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);
    
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
    
  private:
    void initialize();

    double get_model_fraction(int model_size);

    enum ConfigList{
      NoMatch,
      Threshold,
      ModelThreshCalculated
    };
    
    float accThreshold;

    map<string, ConfigList> ConfigMap;
    
    Methods::MDR mdr_calculator;
    
};
}
#endif
