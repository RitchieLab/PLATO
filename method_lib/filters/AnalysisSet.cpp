//AnalysisSet.cpp

#include "AnalysisSet.h"
#include "OutputWriterFactory.h"
using namespace Methods;

namespace Filters{

///
/// Default constructor
///
AnalysisSet::AnalysisSet(){
  outputStream = NULL;
  reducedset = NULL;
}

///
/// Uses job ID information in reporting results
/// @param config ConfigReader with configuration parameters for this set
/// @param outfile string with output file name
/// @param job_id job ID for this set in the database
///
AnalysisSet::AnalysisSet(ConfigReader * config, string outfile){
  outputStream = NULL;
  reducedset = NULL;
  set_parameters(config, outfile);
}


//Destructor
AnalysisSet::~AnalysisSet(){
  for(unsigned int currAnalysis=0; currAnalysis<analyses.size(); currAnalysis++){
    delete analyses[currAnalysis];
  }
  if(outputStream != NULL){
    delete outputStream;
  }
  if(reducedset != NULL){
    delete reducedset;
  }
}



///
/// Sets the dataset for running plato-type analysis
/// @param set DataSet 
/// 
void AnalysisSet::set_dataset(Methods::DataSet* set){
  dataset = set;
  create_reduced_set(dataset);
}



///
/// Sets parameters for run 
/// @param config ConfigReader with the options for the run 
/// @param outfile string with base name for output to file
/// @returns
///
void AnalysisSet::set_parameters(ConfigReader * config, string outfile){

  srand(config->get_seed());

//   set_datafile(config);

  // establish output type 
  outputStream = OutputWriterFactory::create_writer(config->get_output()); 
  if(!config->output_name_set())
    outputStream->establish_stream(config->get_dataset_info(), outfile); 
  else
    outputStream->establish_stream(config->get_dataset_info(), config->get_output_name());

  unsigned int currAnalysis;
  unsigned int numAnalyses = config->num_sets();
  vector<FilterParams>::iterator paramIter;
  vector<FilterParams> params;
  string filterID;
    
  for(currAnalysis=0; currAnalysis < numAnalyses; currAnalysis++){
    params = config->get_params(currAnalysis);
    analyses.push_back(new Analysis);
    for(paramIter = params.begin(); paramIter != params.end();
      paramIter++){
        Filter * newFilter = FilterFactory::create_filter(paramIter->filter_name);
        create_filter(*paramIter, newFilter);
        analyses.back()->add_filter(newFilter); 
        filterID = outputStream->record_parameters(*paramIter);
        newFilter->set_filterID(filterID);
    }
    analyses.back()->add_loci_options(config->get_loci_options(currAnalysis));
  }
}


/// 
/// Sets parameters in filter and creates any subfilters
/// for the filter
/// @param param FilterParams containing all the parameters
/// @param currFilter Filter pointer 
/// @return
/// @throws PlatoExcept on error
///
void AnalysisSet::create_filter(FilterParams & param, Filter * currFilter){
  currFilter->set_params(param.filter_params, dataset);  
  vector<FilterParams>::iterator paramIter;
  for(paramIter = param.sub_filters.begin(); paramIter != param.sub_filters.end();
    paramIter++){
    Filter * subfilter = FilterFactory::create_filter(paramIter->filter_name);
    create_filter(*paramIter, subfilter); 
    currFilter->add_subfilter(subfilter);
  }
} 


///
/// overload << for output
/// @param  os output stream
/// @param  anSet AnalysisSet to output
/// @return  output stream
///
std::ostream & operator << (ostream & os, AnalysisSet & anSet){
  for(unsigned anIndex=0; anIndex < anSet.analyses.size(); anIndex++){
    os << "Analysis #" << anIndex + 1 << endl;
    anSet.analyses[anIndex]->output_analysis(os, *anSet.dataset);
    os << endl;
  }
  return os;
}


///
/// Create a duplicate dataset that only contains valid individuals
/// and markers.  Exclude any marker or sample that is not enabled
/// @param orig_set DataSet
///
void AnalysisSet::create_reduced_set(DataSet* orig_set){
  if(reducedset != NULL)
    delete reducedset;
  reducedset = new DataSet;
  *reducedset = *orig_set;
  
  // clear the vectors
  reducedset->get_markers()->clear();
  reducedset->get_samples()->clear();
  reducedset->get_families()->clear();
  
  // for markers need to use setLoc to get the correct position
  // in the marker_map
  // need to replace this information at the end of the run
  original_loc_indexes.clear();
  
  vector<Marker*>* reduced_markers = reducedset->get_markers();
  vector<Marker*>* orig_markers = orig_set->get_markers();

  
  // for this filtering location is irrelevant so put marker
  // map in same order as markers
  int num_markers = int(orig_set->get_markers()->size());
  for(int i = 0; i<num_markers; i++){
  
    if(orig_markers->at(i)->isEnabled()){
      reduced_markers->push_back(orig_markers->at(i));
//       original_loc_indexes.push_back(orig_markers->at(i)->getLoc());
//       reduced_markers->back()->setLoc(i);
//       reduced_marker_map->push_back(i);
    }
  }
  
  // can just transfer samples
  vector<Sample*>* orig_samples = orig_set->get_samples();
  vector<Sample*>* reduced_samples = reducedset->get_samples();
  
  int num_inds = int(orig_samples->size());
  for(int i=0; i<num_inds; i++){
    if(orig_samples->at(i)->isEnabled()){
      reduced_samples->push_back(orig_samples->at(i));
    }
  }
  
  /// Have to establish the families for new one ///
  reducedset->recreate_family_vector();
  
  // at end need to set affected and unaffected vectors 
  reducedset->set_affection_vectors();
  
  // create marker name map
  reducedset->create_marker_name_map();
  
}


///
/// Returns original set to starting condition.  (use setLoc)
///
void AnalysisSet::return_original_set(){
  int num_markers = int(original_loc_indexes.size());
  for(int i=0; i<num_markers; i++){
    dataset->get_markers()->at(i)->setLoc(original_loc_indexes[i]);
  }
}


///
/// Calculates time estimate
/// @return time in seconds
///
double AnalysisSet::run_time_estimate(){
  double total_time = 0.0;

  for(size_t currAnalysis=0; currAnalysis < analyses.size();
    currAnalysis++){  
      total_time += analyses[currAnalysis]->run_time_estimate(*reducedset);
  }
  return total_time;
}


///
/// runs analysis for this set on each Analysis object contained
/// Outputs results to location specified
/// @return none
/// 
void AnalysisSet::run(){

  for(unsigned int currAnalysis=0; currAnalysis < analyses.size();
    currAnalysis++){
        analyses[currAnalysis]->analyze(*reducedset);
        if(outputStream != NULL)
          analyses[currAnalysis]->output_results(outputStream, *reducedset);
        // this is added for integration in the WASP framework
        // pass original dataset for enabling / disabling loci
        analyses[currAnalysis]->disable_loci(*dataset, *reducedset);
  }
}
}
