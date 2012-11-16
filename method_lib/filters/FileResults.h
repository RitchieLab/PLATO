// FileResults.h

#ifndef __FILERESULTS_H__
#define __FILERESULTS_H__

#include "OutputWriter.h"
#include <fstream>
#include <iomanip>

namespace Filters{
///
/// Implements OutputWriter interface for writing to
/// text file
///


/// Output results to a file
class FileResults: public OutputWriter{
  
  public:
    
    /// Constructor
    FileResults(){};

    /// Closes file
    virtual ~FileResults();
  
    /// Opens output file
    virtual void establish_stream(string outputInfo, string outfile);
    
    /// Adds analysis name as header for adding analysis results
    /// one per line
    virtual void set_analysis(string analysisName);
    
    /// Adds analysis names to output file for display in 
    /// table format.
    virtual void set_analysis(vector<string> analysisNames);
    
    /// Records all results in the specified set as an additional
    /// line in the output file per result
    virtual void record_results(ResultSet & resultset, Methods::DataSet & dataset);
    
    /// Records last result in result set as an additional
    /// line in result file
    virtual void record_last_result(Result & res, string analysisName);
    
    /// Records last result in result set -- assumes that the 
    /// analysis type has been set using the set analysis
    virtual void record_last_result(Result & res);

    /// Currently just returns without putting all
    /// parameters in file
    virtual string record_parameters(FilterParams & params){return "no ID";}
    
    /// Appends results to the file
    void append_results(ResultSet & resultset, Methods::DataSet & dataset);
    
    /// Adds title to first row of output
    void set_title(string title){outfile << "Loci\t" << title << endl;}
  
  private:
    ofstream outfile;
    string current_analysis, outfilename;
    vector<string> analysis_names;
    
};
}

#endif

