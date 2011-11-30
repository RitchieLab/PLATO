// OutputWriter.h

#ifndef __OUTPUTWRITER_H__
#define __OUTPUTWRITER_H__

#include "PlatoDefinitions.h"
#include "Stringmanip.h"
#include <string>
#include "ResultSet.h"
#include <DataSet.h>
#include "FilterDefinitions.h"
#include "PlatoExcept.h"
using namespace std;

///
/// Virtual base class for directing
/// the output from analysis.
///
namespace Filters{
/// Base class for outputting results
class OutputWriter{

  public:
  
  
    /// Virtual destructor
    virtual ~OutputWriter(){}
  
    /// 
    /// Performs basic operations to enable stream
    /// @param outputInfo specifies filename or other information
    /// @param outFile string output filename
    /// @param jobID unsigned int with job ID from PLATO database
    /// @return 
    /// @throws PlatoExcept on error
    ///
    virtual void establish_stream(string outputInfo, string outFile)=0;
    
    ///
    /// Sets analysis type to display as output or for
    /// inserting in database
    /// @param analysisName string containing text name of filter
    /// @return 
    ///
    virtual void set_analysis(string analysisName)=0;
    
    ///
    /// Sets analysis type to display as output or for
    /// inserting in database
    /// @param analysisNames vector of string containing text name of filter
    /// @return 
    ///
    virtual void set_analysis(vector<string> analysisNames)=0;
    
    ///
    /// Records all results in the specified set
    /// @param resultset ResultSet containing all the results
    /// @param dataset PlatoDataset containing loci information
    /// @return 
    /// @throws PlatoExcept on error
    ///
    virtual void record_results(ResultSet & resultset, Methods::DataSet & dataset)=0;
    
    ///
    /// Records last result in result set
    /// @param res Result containing information to record
    /// @param analysisName string giving name of analysis
    /// @return 
    /// @throws PlatoExcept on error
    /// 
    virtual void record_last_result(Result & res, string analysisName)=0;
    
    ///
    /// Records last result in result set -- assumes that the 
    /// analysis type has been set using the set analysis
    /// function
    /// @param res Result containing information to record
    /// @return 
    /// @throws PlatoExcept on error
    /// 
    virtual void record_last_result(Result & res)=0;

    ///
    /// Records parameters for the specified analysis
    /// @param params PARAMS for the analysis
    /// @return filter ID 
    /// @throws PlatoExcept on error
    ///
    virtual string record_parameters(FilterParams & params)=0;
  
  private:
  
};
}
#endif
