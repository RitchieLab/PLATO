//MarkerGenoFilter.h

#ifndef __MARKERGENOFILTER_H__
#define __MARKERGENOFILTER_H__

#include "Filter.h"
#include <MarkerGenoEff.h>
#include <time.h>

///
/// Filter conducts marker efficiency calculation.
///

namespace Filters{

/// Runs Marker genotyping efficiency
class MarkerGenoFilter: public Filter{
  
  public:
    MarkerGenoFilter();
    MarkerGenoFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);

    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);

  private:
    void initialize();

    enum ConfigList{ 
      NoMatch,
      Threshold
    };

    float efficiencyThreshold;
    Methods::MarkerGenoEff eff_calculator;
    map<string, ConfigList> ConfigMap;

};
}
#endif
