//Analysis.h

#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <vector>
#include "Filter.h"
#include <iomanip>
#include <filters/ResultsManager.h>

namespace Filters{
/// Contains filters for running analysis
class Analysis{
  
  public:
    Analysis();
    Analysis(const Analysis & origAnalysis);
    virtual ~Analysis();
    Analysis & operator=(const Analysis & origAnalysis);
    
    friend ostream & operator<<(std::ostream & os, Analysis & an);
    
    void output_results(OutputWriter * resultStream, Methods::DataSet & dataset);
    
    void output_analysis(std::ostream & os, 
      Methods::DataSet & dataset);

    virtual void analyze(Methods::DataSet & dataset);
    void add_filter(Filter * newFilter);
//     void add_filter(std::vector<string> & filterNames, 
//       std::vector<std::map<string, string> > & params);  
    void add_loci_options(vector<ListOptions> lociOptions){locSizeOptions = lociOptions;}
  
    // Disables missing loci in dataset for passing back to wasp 
    void disable_loci(Methods::DataSet& orig_dataset, Methods::DataSet& analysis_dataset);
  
    void set_total_proc(int num){total_proc = num;}
    int get_total_proc(){return total_proc;}
    void set_proc_num(int num){proc_num = num;}
    int get_proc_num(){return proc_num;}
    
    double run_time_estimate(Methods::DataSet& dataset);
  
//  protected:
  
    void initialize();
  
    std::vector<Filter *> filters;
    std::vector<ListOptions> locSizeOptions;
    ResultsManager * resultHandler;
  
    int proc_num, total_proc;
  
};
}
#endif
