//MDRPDTFilter.h

#ifndef __MDRPDTFILTER_H__
#define __MDRPDTFILTER_H__

#include "Filter.h"
#include <MDRPDT.h>

///
/// Filter conducts MDR PDT
///
namespace Filters{
/// Runs MDR PDT
class MDRPDTFilter: public Filter{
  
  public:
    MDRPDTFilter();
    MDRPDTFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);

    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
    
  private:
    void initialize();

    enum ConfigList{
      NoMatch,
      Threshold,
      RandSeed,
      Ptests,
      Crossvals,
      MinComb,
      MaxComb
    };
    
    float threshold;

    map<string, ConfigList> ConfigMap;
    
    int n_perms;
    
    Methods::MDRPDT mdrpdt_calculator;
    
};
}
#endif
