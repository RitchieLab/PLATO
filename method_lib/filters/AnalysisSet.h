//PlatoConfig.h

#ifndef __ANALYSISSET_H__
#define __ANALYSISSET_H__

#include <vector>
#include <iostream>
#include <DataSet.h>
#include "Analysis.h"
#include <fstream>
#include "ConfigReader.h"

///
/// Contains all analyses in the set <br>
/// An Analysis object is needed if the analyses
/// are to be run in parallel
///
namespace Filters{

/// Groups analysis objects for conducting analysis on same dataset
class AnalysisSet{
  
  public:
    AnalysisSet();
    AnalysisSet(ConfigReader * config, string outfile);
    virtual ~AnalysisSet();
    
    virtual void run();
    //void read_config(string configFile);
    friend std::ostream & operator << (std::ostream & os, AnalysisSet & anset);


    /// sets parameters for run
    virtual void set_parameters(ConfigReader * config, string outfile);
    
    /// sets the dataset for running
    void set_dataset(Methods::DataSet* set);
    
    /// gets time estimate
    double run_time_estimate();

  protected:  
    /// Set filters parameters and create subfilters  
    void create_filter(FilterParams & param, Filter * currFilter);
    
    /// Create a reduced dataset to run through plato
    void create_reduced_set(Methods::DataSet* orig_set);
    
    /// Returns original set to condition (use setLoc)
    void return_original_set();
    
    Methods::DataSet* dataset, *reducedset;
    
    OutputWriter * outputStream;
    std::vector<Analysis *> analyses;
    std::vector<int> original_loc_indexes;
};

}
#endif
