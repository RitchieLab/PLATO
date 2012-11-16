//ContingencyFilter.h

#ifndef __CONTINGENCYFILTER_H__
#define __CONTINGENCYFILTER_H__

#include "Filter.h"
#include "Permutation.h"
#include <map>

namespace Filters{
/// Filter calculates stats on contingency tables
class ContingencyFilter: public Filter{

  public:
    ContingencyFilter();
    ContingencyFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);   

    /// returns estimate of run time for this filter
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
    
  private:
    void initialize();
  
    void set_text_outname();
  
    // add to list below for each additional parameter
    enum ConfigList{ 
      NoMatch,
      Threshold,
      Permutations,
      ReportLoci
    };
    
    bool perm_completed;
    map<unsigned int, bool> report_loci;
    float threshold;
    unsigned int n_perms;
    Permutation perm;
    map<string, ConfigList> ConfigMap;

};

}
#endif
