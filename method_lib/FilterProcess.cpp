//FilterProcess.cpp

#include "FilterProcess.h"
#include "filters/ConfigFileReader.h"
#include "filters/AnalysisSet.h"

using namespace Filters;

namespace Methods{
///
/// Basic constructor
///
FilterProcess::FilterProcess(){
  dataset = NULL;
  initialize();
}

///
/// Sets dataset and then initializes
///
FilterProcess::FilterProcess(DataSet* ds){
  resetDataSet(ds);
  initialize();
}

///
/// Sets configuration filename, dataset and then initializes
///
FilterProcess::FilterProcess(DataSet* ds, string ConfigFileName){
  resetDataSet(ds);
  setConfigFilename(ConfigFileName);
  initialize();
}


///
/// Sets initial state
///
void FilterProcess::initialize(){
  configfilename = "";
  outputname = "plato_filter.out";
}


///    
/// Returns estimated runtime in seconds
/// @return estimated run time in seconds
///
double FilterProcess::GetEstimatedRunTime(DataSet* set, string ConfigFileName){
  resetDataSet(set);
  setConfigFilename(ConfigFileName);
  return GetEstimatedRunTime();
}

///    
/// Returns estimated runtime in seconds
/// @return estimated run time in seconds
///
double FilterProcess::GetEstimatedRunTime(){
  AnalysisSet analysis_set;
  analysis_set.set_dataset(dataset);
  ConfigFileReader config;
  config.read_config(configfilename);
  
  analysis_set.set_parameters(&config, outputname);
  
  return analysis_set.run_time_estimate();
}
 
///
/// Runs process -- at end any loci that are excluded by filtering will be disabled
/// @param set DataSet
/// @param string
///
void FilterProcess::calculate(DataSet* set, string ConfigFileName){
  resetDataSet(set);
  setConfigFilename(ConfigFileName);
  calculate();
}


///
/// Runs process -- at end any loci that are excluded by filtering will be disabled.
/// Runs on the DataSet and configuration file name already set.
///
void FilterProcess::calculate(){
  AnalysisSet analysis_set;
  analysis_set.set_dataset(dataset);
  ConfigFileReader config;
  config.read_config(configfilename);
  analysis_set.set_parameters(&config, outputname);
  analysis_set.run();
}

///
/// sets DataSet 
/// @param 
///
void FilterProcess::resetDataSet(DataSet* ds){
  dataset = ds;
}
}

