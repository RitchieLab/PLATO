//FilterProcess.h

#ifndef __FILTERPROCESS_H__
#define __FILTERPROCESS_H__

#include "DataSet.h"
#include "StepOptions.h"
#include "MethodException.h"

///
/// Method that replaces the standard standalone PLATO process.
///

/// Method class that performs MDR
namespace Methods{
class FilterProcess{

  public:
    FilterProcess();
    FilterProcess(DataSet* ds);
    FilterProcess(DataSet* ds, string ConfigFileName);
    
    /// Returns estimated runtime in seconds
    double GetEstimatedRunTime(DataSet* set, string filename);
    
    /// Returns estimated runtime in seconds
    double GetEstimatedRunTime();
    
    /// Runs process -- at end any loci that are excluded by filtering will be disabled
    void calculate(DataSet* set, string filename);
    
    /// Runs process on DataSet and configuration filename already set in class
    void calculate();
    
    /// Sets configuration file name
    inline void setConfigFilename(string filename){configfilename = filename;}
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// 
    inline void setOutputName(string out_name){outputname = out_name;}

  private:
  
    void initialize();
    
    DataSet* dataset;
  
    string configfilename, outputname;
    

};
}
#endif

