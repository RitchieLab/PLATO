//FileResults.cpp
#include "FileResults.h"

using namespace Methods;
namespace Filters{
///
/// Destructor closes output file if open
/// 
FileResults::~FileResults(){
  if(outfile.is_open()){
    outfile.flush();
    outfile.close();
  }
}


/// 
/// Opens output file
/// @param outputFileName specifies filename or other information
/// @param outFile string output name
/// @param jobID ID number from PLATO database for the present job
/// @return
/// @throws PlatoExcept on error
///
void FileResults::establish_stream(string outputFileName, string outFile){
  outputFileName = outFile;
  outfile.open(outputFileName.c_str(), ios::out);

  if(!outfile.is_open()){
    throw PlatoExcept(outputFileName + " unable to open for writing for results");
  }
//   outfile << "Database Job ID=" << jobID << endl;
  outfile << "Filter process output" << endl;
  outfilename = outFile;
}


///
/// Sets analysis info and outputs name as header
/// @param analysisName string containing text name of filter
/// @return 
///
void FileResults::set_analysis(string analysisName){
  current_analysis = analysisName;
  outfile << "Analysis: " << current_analysis << "\n";
}


///
/// Sets analysis info and outputs name as header
/// @param analysisNames vector containing text names of analyses 
///        used and should be in order that results appear
///        in the result set
/// @return 
///
void FileResults::set_analysis(vector<string> analysisNames){
  analysis_names = analysisNames;
  // set up header row for output format
  outfile  << "Loci";
  for(unsigned int aname=0; aname < analysis_names.size(); aname++){
    outfile << "\t" << analysis_names[aname];
  }
  outfile << endl;
}


///
/// Records last result in result set as an additional
/// line in result file -- used when needing to record
/// all results during the analysis so that even results
/// that are later trimmed are captured
/// @param res Result containing information to record
/// @param analysisName string giving name of analysis
/// @return 
/// @throws PlatoExcept on error
/// 
void FileResults::record_last_result(Result & res, string analysisName){
  record_last_result(res);
}


///
/// Records last result in result set -- assumes that the 
/// analysis type has been set using the set analysis
/// function -- used when needing to record
/// all results during the analysis so that even results
/// that are later trimmed are captured
/// @param res Result containing information to record
/// @ret none
/// @throws PlatoExcept on error
/// 
void FileResults::record_last_result(Result & res){
  if(res.analysisScores.back() >=0)
    outfile << setw(10) << res.analysisScores.back() << "\n";
  else
    outfile << setw(10) << " " << "\n";
}


///
/// Appends results to the file
/// @param resultset ResultSet containing all the results
/// @param dataset DataSet containing info on the loci
///
void FileResults::append_results(ResultSet & resultset, Methods::DataSet & dataset){
  if(outfile.is_open()){
    outfile.close();
  }
  
  outfile.open(outfilename.c_str(), ios::app);
  record_results(resultset, dataset);
  outfile.close();
  
}


///
/// Records all results in the specified set as an additional
/// line in the output file per result
/// @param resultset ResultSet containing all the results
/// @param dataset DataSet containing info on the loci
/// @return 
/// @throws PlatoExcept on error
///
void FileResults::record_results(ResultSet & resultset,
  DataSet & dataset){
  
  ResultIter resIter;
  unsigned int currLoc, currValue, totalValues;
  string loci;
  std::_Ios_Fmtflags oldflags = outfile.setf(ios::left);
  
  // determine number of values in each result
  // all will have same number
  resIter = resultset.begin();
  totalValues = resIter->analysisScores.size();
  vector<Marker*> markers = *(dataset.get_markers());
  vector<int> marker_map = *(dataset.get_marker_map());
  
  for(; resIter != resultset.end(); resIter++){
    // output loci combination for result
    loci = "";
    for(currLoc = 0; currLoc < resIter->genoCombination.size(); currLoc++){
//       loci += dataset.get_locus_name(resIter->genoCombination[currLoc]) + " ";
//       loci += markers[marker_map[resIter->genoCombination[currLoc]]]->getRSID() +
//         " ";
      loci += markers[resIter->genoCombination[currLoc]]->getRSID() + " ";
    }
    outfile << setw(10) << loci;
    for(currValue = 0; currValue < totalValues; currValue++){
      if(resIter->analysisScores[currValue] > -1e30)
        outfile << "\t" << resIter->analysisScores[currValue];
      else
        outfile << "\t";
    } 
    outfile << endl;
  }
  
  outfile.setf(oldflags);
}
}


