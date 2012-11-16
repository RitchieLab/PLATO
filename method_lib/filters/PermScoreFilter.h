//PermScoreFilter.h

#ifndef __PERMSCOREFILTER_H__
#define __PERMSCOREFILTER_H__

///
/// Runs sub filters and creates a distribution for
/// each snp and each sub filter that can be used for
/// assigning p values.
///

#include "Filter.h"
#include "Permutation.h"
#include <map>
#include <set>
namespace Filters{
/// Filter calculates stats on contingency tables
class PermScoreFilter: public Filter{

  public:
    PermScoreFilter();
    PermScoreFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);
    
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
    
    vector<vector<multiset<float> > >distributions;
    
};
}



#endif
