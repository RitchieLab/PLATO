//MaxScoreFilter.h

#ifndef __MAXSCOREFILTER_H__
#define __MAXSCOREFILTER_H__

///
/// Runs sub filters and then selects single best score from those
/// (based on p value).  Can then permute data and create p value 
/// score based on that distribution
///

#include "Filter.h"
#include "Permutation.h"
#include <map>
#include <set>
namespace Filters{
/// Filter calculates stats on contingency tables
class MaxScoreFilter: public Filter{

  public:
    MaxScoreFilter();
    MaxScoreFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);
 
    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);
 
  private:
    void initialize();
  
    void run_perms(Methods::DataSet & data);
  
    void set_text_outname();
  
    // add to list below for each additional parameter
    enum ConfigList{ 
      NoMatch,
      Threshold,
      Permutations
    };
    
    bool perm_completed, name_adjusted;
    unsigned int n_perms;
    float threshold;
    Permutation perm;
    
//     map<unsigned int, bool> report_loci;
//     float threshold;
    
    map<string, ConfigList> ConfigMap;
    vector<multiset<float> > distributions;
    
};

}


#endif
