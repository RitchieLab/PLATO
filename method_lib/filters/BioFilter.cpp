//BioFilter.cpp

#include "BioFilter.h"
#include "FileResults.h"
#include "BioBinaryReader.h"

using namespace Methods;
namespace Filters{
///
/// constructor -- initialize variables
///
BioFilter::BioFilter():Filter("BioFilter"){
  initialize();
}


///
/// Alternative constructor that sets filter name
/// @param filterName name of filter
///
BioFilter::BioFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void BioFilter::initialize(){
  ConfigMap["TABLE"] = TableName;
  ConfigMap["INFILE"] = InFileName;
  ConfigMap["TABLEOUT"] = TableOut;
  ConfigMap["INTERACTIONINCLUDED"] = InteractionIncluded;
  ConfigMap["COVARINCLUDED"] = CovariateIncluded;
  ConfigMap["OUTSCHEMA"] = OutSchema;
  ConfigMap["INSCHEMA"] = InSchema;
  ConfigMap["MAXMODELSIZE"] = MaxModelSize;
  ConfigMap["OUTBIOFILE"] = OutTextFile;
  interact_included = true;
  include_cov = false;
  in_schema = "Analysis";
  out_schema = "Analysis";
  max_model_size = 2;
  in_filename = "";
  out_text_file = "";
}


///
/// Fills unsigned int vector with indexes taken from 
/// string ids
/// @param index_combo
/// @param string_combo
/// @param dataset
///
void BioFilter::fill_index(vector<unsigned int>& index_combo, vector<string>& string_combo,
  DataSet& dataset){
  
  for(unsigned int i=0; i < string_combo.size(); i++){
    index_combo[i] = dataset.get_locus_index(string_combo[i]);
  }

}


string BioFilter::current_time(){
time_t rawtime;
struct tm * timeinfo;
time ( &rawtime );
timeinfo = localtime ( &rawtime ); 
 return asctime(timeinfo);
}


///
/// Calculates contingency tables for all
/// subfilters that use the contingency table class.
/// Optimized for quicker speed.  Currently can
/// only be used with the text output.
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void BioFilter::analyze(ResultSet & resultList, DataSet & dataset){

  if(resultList.size() > 0){
    analyze_list(resultList, dataset);
  } 
#if MYSQLDEFINED
  else if(bio_tables.size() > 0){
    analyze_from_db(resultList, dataset);
  }
#endif
  else if(in_filename.length() > 0){
    analyze_from_text(resultList, dataset);
  }
  
}


///
/// Run list of combinations acquired from binary file produced
/// by Eric's app.
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
///
void BioFilter::analyze_from_text(ResultSet & resultList, DataSet & dataset){

  unsigned int curr_sub;
  unsigned int num_subs = subfilters.size();

  set_text_outname();
  
  vector<vector<string> > snp_ids;
  bool combos_complete;
  unsigned int curr_combo;
  vector<unsigned int> index_combo;
  
  Result newResult;
  ResultSet temp_list;  

  BioBinaryReader binreader(in_filename);
  FileResults recorder;

  recorder.establish_stream(out_text_file, out_text_file);

  recorder.set_title(name);
  
  do{
    temp_list.clear();

    combos_complete = binreader.get_combinations(snp_ids);
    
    if(snp_ids.size() > 0){
      for(curr_combo = 0; curr_combo < snp_ids.size(); curr_combo++){
 
 
        Result newResult;
        newResult.genoCombination.resize(snp_ids[curr_combo].size());
 
        // if catch exception indicating missing a locus name
        // continue with next combination
        try{
          fill_index(newResult.genoCombination, snp_ids[curr_combo], dataset);
        }catch(MethodException& de){
          cout << "Model skipped: ";
          for(unsigned int z=0; z<snp_ids[curr_combo].size(); z++){
            cout << snp_ids[curr_combo][z] << " ";
          }
          cout << endl;
          continue;
        }
        temp_list.push_back(newResult);
      }
    }   
    for(curr_sub=0; curr_sub < num_subs; curr_sub++){
      subfilters[curr_sub]->analyze(temp_list, dataset);
    }
       
    // add to end of actual list
    resultList.splice(resultList.end(), temp_list);
    temp_list.clear();
 
    // insert any results
    recorder.append_results(resultList, dataset);
    
    // clear list so that list will not be retained
    resultList.clear();
      
    // clear IDs and try to get more
    snp_ids.clear();
  
  }while(!combos_complete);
  
}


///
/// Fills result list with 100 models for running sub filters
/// and determining time for run
/// @param testList 
/// @param dataset DataSet
///
void BioFilter::fill_list(ResultSet& templist, Methods::DataSet& dataset){
  int total_loci = dataset.num_loci();
  int largest_size = 2;
  
  vector<int> selected;
  int randloc;
  while(selected.size() < 10){
    randloc = int(float(rand()) / RAND_MAX * total_loci);
    if((dataset.get_markers())->at(randloc)->isEnabled()){
      selected.push_back(randloc);
    }
  }
  Result newResult;
  newResult.genoCombination.resize(largest_size);
  
  // create based on largest size of model to be checked
  int startloc = 0;
  size_t currloc;
  for(startloc=0; startloc < 10; startloc++){
    currloc = startloc;
    Result newResult;
    for(int i=0; i<largest_size; i++){
      newResult.genoCombination.push_back(selected[currloc]);
      ++currloc;
      if(currloc > selected.size()-1){
        currloc = 0;
      }
    }
    templist.push_back(newResult);
  }
  
 
  ResultSet addlist;
  addlist = templist;
//  templist.push_back(templist);
  for(int i=0; i<9; i++)
    templist.insert(templist.end(), addlist.begin(), addlist.end());
}


///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate BioFilter::estimate_run_time(double num_models, Methods::DataSet& dataset,
  ResultSet& testList, int model_size){
 
  // how about using a binreader to read number of models and fill initial list
  BioBinaryReader binreader(in_filename);
  unsigned int nmodels=0;
  std::vector<std::vector<std::string> > combos;
  bool done=false;
  do{
    done = binreader.get_combinations(combos);
    nmodels += combos.size();
    combos.clear();
  }while(!done);

  // if result list size is empty insert 100 models into list for running
  // sub filters
  if(testList.size() == 0)
    fill_list(testList, dataset);

  Filter::ProcessEstimate time_estimate;
  num_models = nmodels;

  Filter::ProcessEstimate sub_estimate;
  sub_estimate.num_models = num_models;
  double time_used = 0.0;
  
  // need to call estimate on each sub filter
  for(size_t curr_sub=0; curr_sub < subfilters.size(); curr_sub++){
    ResultSet subList = testList;
    sub_estimate = subfilters[curr_sub]->estimate_run_time(sub_estimate.num_models, 
      dataset, subList, 2);
    time_used += sub_estimate.seconds;
  }

  time_estimate.seconds = time_used;
  // all models cleared at end
  time_estimate.num_models = 0;
  
  return time_estimate;
}



///
/// Run list of combinations passed to filter
/// and store in database all the results
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void BioFilter::analyze_list(ResultSet & resultList, DataSet & dataset){
  
  unsigned int curr_sub;
  unsigned int num_subs = subfilters.size();
  
  // run sub filters on list
  for(curr_sub=0; curr_sub < num_subs; curr_sub++){
    subfilters[curr_sub]->analyze(resultList, dataset);
  }
  
  
#ifdef MYSQLDEFINED  
  // output any remaining results to the database table
  // set connection parameters for these results
  recorder.set_connect_params("munster.mc.vanderbilt.edu", "plato", "genomewide", out_schema);
  
  // check for presence of output table and create if necessary
  recorder.check_table(max_model_size, interact_included, output_table_name, out_schema);
  
  string output_table_full = out_schema + "." + output_table_name;
  
  recorder.record_results(resultList, dataset, output_table_name, interact_included, include_cov);
#endif

  // clear list so that list will not be retained
  resultList.clear();
  
}


#ifdef MYSQLDEFINED
///
/// Get combinations from the database or binary file and do analysis
/// @param  resultList ResultSet to store results
/// @param  dataset DataSet to process
/// @return  none
///
void BioFilter::analyze_from_db(ResultSet & resultList, DataSet & dataset){

  unsigned int curr_sub;
  unsigned int num_subs = subfilters.size();

  set_text_outname();

  vector<table_info>::iterator table_iter;
  
  vector<vector<string> > snp_ids;
  bool combos_complete;
  unsigned int curr_combo;
  vector<unsigned int> index_combo;
  
  Result newResult;
  ResultSet temp_list;

  recorder.set_connect_params("munster.mc.vanderbilt.edu", "plato", "genomewide", out_schema);

  // check for presence of output table and create if necessary
  recorder.check_table(max_model_size, interact_included, output_table_name, out_schema);
  
  string output_table_full = out_schema + "." + output_table_name;

  BioFilterCombinations bio_combo("munster.mc.vanderbilt.edu", "plato", "genomewide", in_schema);

  // do one table/view at a time
  for(table_iter = bio_tables.begin(); table_iter != bio_tables.end(); table_iter++){
    bio_combo.fill_combinations(table_iter->name, table_iter->combo_size);
    combos_complete = bio_combo.get_combinations(snp_ids);
    
    newResult.genoCombination.resize(table_iter->combo_size);

    while(snp_ids.size() > 0){
      
      // fill result list
      for(curr_combo = 0; curr_combo < snp_ids.size(); curr_combo++){
        // if catch exception indicating missing a locus name
        // continue with next combination
        try{
          fill_index(newResult.genoCombination, snp_ids[curr_combo], dataset);
        }catch(MethodException& de){
          continue;
        }
        temp_list.push_back(newResult);
      }
      // run sub filters on list
      for(curr_sub=0; curr_sub < num_subs; curr_sub++){
        subfilters[curr_sub]->analyze(temp_list, dataset);
      }

      // add to end of actual list
      resultList.splice(resultList.end(), temp_list);
      temp_list.clear();
 
      // insert any results
      recorder.record_results(resultList, dataset, output_table_full, interact_included, include_cov);
      // clear list so that list will not be retained
      resultList.clear();
      
      // clear IDs and try to get more
      snp_ids.clear();
      if(!combos_complete)
        combos_complete = bio_combo.get_combinations(snp_ids);
    }
  } 

}

#endif

///
/// Sets name of this filter to reflect sub filters
/// used -- Will appear as column headers in text file
/// @return
///
void BioFilter::set_text_outname(){
  if(subfilters.size() == 0)
    return;
  name = subfilters[0]->getName();
  unsigned int num_subs = subfilters.size();
  for(unsigned int curr_sub=1; curr_sub < num_subs; curr_sub++){
      name += "\t" + subfilters[curr_sub]->getName();
  }  
}


///
/// Sets parameters for the filter and throws a FilterExcept when 
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @param  set DataSet
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void BioFilter::set_params(PARAMS params, DataSet* set){

  std::map<string, string>::iterator configIter;
  stringstream ss;
  table_info table;
  string value;
  
  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){
    switch(ConfigMap[configIter->first]){
      case TableName:
        
        ss.str(configIter->second);
        while(!ss.eof()){
          ss >> table.combo_size;
          ss >> table.name;
          bio_tables.push_back(table);
        }
        break;
      case TableOut:
        output_table_name = configIter->second;
       break;
      case InSchema:
        in_schema = configIter->second;
        break;
      case OutSchema:
        out_schema = configIter->second;
        break;
      case InFileName:
        in_filename = configIter->second;
        break;
      case OutTextFile:
        out_text_file = configIter->second;
        break;
      case InteractionIncluded:
        value = Stringmanip::to_upper(configIter->second);
        if(value.compare("TRUE") == 0)
          interact_included = true;
        else
          interact_included = false;
        break;
      case CovariateIncluded:
         value = Stringmanip::to_upper(configIter->second);
        if(value.compare("TRUE") == 0)
          include_cov = true;
        else
          include_cov = false;
        break;       
      case MaxModelSize:
        max_model_size = Stringmanip::stoi(configIter->second);
        break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " + 
          name + " filter\n\n");
    };
  }  
  set->create_marker_name_map();
  
}
}

