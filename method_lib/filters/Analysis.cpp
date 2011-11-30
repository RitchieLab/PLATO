//Analysis.cpp

#include "Analysis.h"
#include <set>

using namespace Methods;
namespace Filters{
///
/// Default Constructor
///
Analysis::Analysis(){
  initialize();
  resultHandler = new ResultsManager;
}

///
/// Copy Constructor for analysis class
/// @param origAnalysis original Analysis object
///
Analysis::Analysis(const Analysis & origAnalysis){
  resultHandler = new ResultsManager;
  *resultHandler = *(origAnalysis.resultHandler);

  for(unsigned int currFilter=0; currFilter < origAnalysis.filters.size();
    currFilter++){
      Filter * newFilter = FilterFactory::create_filter(origAnalysis.filters[currFilter]->getName());
      *newFilter = *origAnalysis.filters[currFilter];
      filters.push_back(newFilter);
  }
}

///
/// operator = overloaded for Analysis object
/// @param origAnalysis source Analysis object
/// @return this object
///
Analysis & Analysis::operator=(const Analysis & origAnalysis){
 
  resultHandler = new ResultsManager;
  *resultHandler = *(origAnalysis.resultHandler);
  
  // clear any existing Filters
  for(unsigned int currFilter=0; currFilter < filters.size(); currFilter++){
    delete filters[currFilter];
  }
  // create copies of filters for this object
  for(unsigned int currFilter=0; currFilter < origAnalysis.filters.size();
    currFilter++){
      Filter * newFilter = FilterFactory::create_filter(origAnalysis.filters[currFilter]->getName());
      *newFilter = *origAnalysis.filters[currFilter];
      filters.push_back(newFilter);
  }
  return *this;
}


///
/// Called by constructors to initialize object
/// Arg: none
/// Ret: none
///
void Analysis::initialize(){
}


// Destructor
Analysis::~Analysis(){
  // destroy all filters in this analysis
  for(unsigned int filterNum=0; filterNum < filters.size(); filterNum++){
    delete filters[filterNum];
  }
}

///
/// Adds filters to the Analysis object
/// @param filterNames names of filters to add
/// @param params maps containing paramters for this analysis
/// @return none
///
// void Analysis::add_filter(std::vector<string> & filterNames, 
//   std::vector<std::map<string, string> > & params){
//     // resize the dynamic array
//     Filter * newFilter;
//     for(unsigned int currFilter=0; currFilter < filterNames.size(); currFilter++){
//       newFilter = FilterFactory::create_filter(filterNames[currFilter]);
//       newFilter->set_params(params[0]);
//       filters.push_back(newFilter);
//     }
// }

///
/// Adds a filter object to vector
/// @param  newFilter new filter in analysis
/// @return  none
///
void Analysis::add_filter(Filter * newFilter){
  filters.push_back(newFilter);
}


///
/// Runs analysis with results stored in list
/// @param dataset for analysis
/// @return none
///
void Analysis::analyze(DataSet & dataset){


  for(unsigned int currFilter=0; currFilter < filters.size(); currFilter++){
    resultHandler->run_filter(filters[currFilter], dataset, locSizeOptions[currFilter]);
  }
}


///
/// Runs time estimate on filters
/// @param dataset DataSet
/// @return time in seconds
///
double Analysis::run_time_estimate(DataSet& dataset){

  Filter::ProcessEstimate estimate;
  double time_estimate = 0.0;
  estimate.num_models = 0;
  estimate.seconds = 0.0;

  for(unsigned int currFilter=0; currFilter < filters.size(); currFilter++){
    estimate = resultHandler->run_filter_time(filters[currFilter], dataset, locSizeOptions[currFilter],
      estimate.num_models);
    time_estimate += estimate.seconds;
  }
  
  return time_estimate;
}




///
/// overload << for output
/// @param  os output stream
/// @param  an Analysis to output
/// @return  output stream
///
std::ostream & operator<<(std::ostream & os, Analysis & an){

  // output all filters names across top
  std::_Ios_Fmtflags oldflags = os.setf(ios::left);
  string loc = "loci";
  os << setw(8) << loc;
  for(unsigned int anIndex=0; anIndex < an.filters.size(); anIndex++){
    os << setw(10) << an.filters[anIndex]->getName();
  }
  os << endl;


  unsigned int currLoc, anIndex;
  
  // output all results
  for(ResultIter currResult = an.resultHandler->resultList.begin(); 
    currResult != an.resultHandler->resultList.end(); 
    ++currResult){
    string loci;
    for(currLoc=0; currLoc < currResult->genoCombination.size();
      currLoc++){
      loci += Stringmanip::itos(currResult->genoCombination[currLoc]) + " ";
    }
    os << setw(8) << loci; 
    
    // output each analysis result
    for(anIndex=0; anIndex < an.filters.size(); anIndex++){
      os << setw(10) << currResult->analysisScores[anIndex];
    }
    
    os << endl; 
  }
  
  os.setf(oldflags);
  return os;
}


///
/// outputs analysis results with locus rs numbers to 
/// supplied output stream
/// @param  os output stream
/// @param dataset DataSet containing rs numbers for the loci
/// @return none
///
void Analysis::output_analysis(std::ostream & os, 
  DataSet & dataset){
  // output all filters names across top
  std::_Ios_Fmtflags oldflags = os.setf(ios::left);
  string loc = "loci";
  os << setw(8) << loc;
  unsigned int anIndex, currLoc;
  for(anIndex=0; anIndex < filters.size(); anIndex++){
    os << setw(10) << filters[anIndex]->getName();
  }
  os << endl;

  string startloc;

  vector<Marker*> markers = *(dataset.get_markers());

  // output all results
  for(ResultIter currResult = resultHandler->resultList.begin(); 
    currResult != resultHandler->resultList.end(); ++currResult){
    string loci;
    for(currLoc=0; currLoc < currResult->genoCombination.size();
      currLoc++){

      loci += markers[currLoc]->getRSID() + " ";
    }
    os << setw(8) << loci; 
    
    // output each analysis result
    for(anIndex=0; anIndex < filters.size(); anIndex++){
      os << setw(10) << currResult->analysisScores[anIndex];
    }
    
    os << endl; 
  }

  os.setf(oldflags);
}



///
/// outputs analysis results to indicated OutputWriter object
/// @param  resultStream OutputWriter specifying location
///         and type of reports
/// @param  dataset PlatoDataset for use in getting loci
///         names
/// @return none
///
void Analysis::output_results(OutputWriter * resultStream, 
  DataSet & dataset){
  vector<string> analysis_names;
  for(unsigned int aIndex=0; aIndex < filters.size(); aIndex++){
    analysis_names.push_back(filters[aIndex]->getName());
  }
  resultStream->set_analysis(analysis_names);
  resultStream->record_results(resultHandler->resultList, dataset);
}



///
/// Sets the loci in the dataset to be enabled if they appear in 
/// result list.  Then clears the result list to free memory.
/// @param datset DataSet
///
void Analysis::disable_loci(DataSet& orig_dataset, DataSet& analysis_dataset){

  // store loci in set to indicate which are present in result list
  // store as string so that no confusion over which are enabled or 
  // not
  set<string> included_loci;
  
  ResultIter resIter;
  unsigned int currLoc;
  vector<Marker*>* markers = orig_dataset.get_markers();
  vector<Marker*>* analysis_markers = analysis_dataset.get_markers();

  for(resIter = resultHandler->resultList.begin(); resIter != resultHandler->resultList.end();
    ++resIter){
    for(currLoc = 0; currLoc < resIter->genoCombination.size(); currLoc++){
      included_loci.insert(analysis_markers->at(resIter->genoCombination[currLoc])->getRSID());
    }
  }

  // iterate through each marker and check to see if it appears in set
  // if not, change enabled to false
  for(unsigned int currMark=0; currMark < markers->size(); currMark++){
    if(included_loci.find(markers->at(currMark)->getRSID()) == included_loci.end()){
      markers->at(currMark)->setEnabled(false);
    }
  }

}

}
