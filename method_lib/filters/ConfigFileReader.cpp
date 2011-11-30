//ConfigFileReader.cpp

#include "ConfigFileReader.h"
#include <sstream>

namespace Filters{

///
/// Alternative constructor that reads in 
/// config file
/// @param config_file string name of configuration file
///
ConfigFileReader::ConfigFileReader(string config_file){
  initialize_keywords();
  read_config(config_file);
}


///
/// Initialize keyword map
/// @return 
///
void ConfigFileReader::initialize_keywords(){
  keywordMap["ANALYSIS"] = keyAnalysis;
  keywordMap["START"] = keyStart;
  keywordMap["END"] = keyEnd;
  keywordMap["DATASET"] = keyDataset;
  keywordMap["SUB"] = keySub;
  keywordMap["RANDSEED"]= keyRandSeed;
  keywordMap["OUT"] = keyOut;
  keywordMap["MISSINGVALUE"] = keyMissingValue;
  keywordMap["LOCICOMBO"] = keyLocCombo;
  keywordMap["INPUT"] = keyDatasetType;
  keywordMap["MAPFILE"] = keyMapFile;
  keywordMap["PSEUDOCONTROLS"] = keyPseudoGenerate;
  keywordMap["OUTNAME"] = keyOutputName;

  lociKeyMap["ADD"] = addLoci;
  lociKeyMap["PASS"] = passLoci;
  lociKeyMap["REMOVE"] = removeLoci;
  
  outputParam = "TEXT";
  dataset_type = "TEXT";
  
  loci_info_id = "";
  output_name = "";
  loci_info_available = false;  
  output_specified = false;
  pseudo_generate = false;
}


///
/// Reads configuration file and stores 
/// parameters for later use
/// @param  configfile name of config file
/// @return  none
/// @throws  PlatoExcept for out of range parameters, missing file, etc.
///
void ConfigFileReader::read_config(string configfile){
  std::ifstream dataStream(configfile.c_str(), ios::in);
  if(!dataStream.is_open()){
    throw PlatoExcept("Failed in attempt to open file " + configfile + "\n");
  }
  
  string line, keyword;
  vector<FilterParams> currAnalysis;

  bool priorAnalysis = false;
  vector<ListOptions> currListOptions;
  ListOptions defaultOptions;
  defaultOptions.op = PassList;
  string pseudo_option;
  
  while(!dataStream.eof()){
    getline(dataStream, line);
   
    if(skip_line(line)){
      continue;
    }
    
    stringstream ss(line);
    ss >> keyword;
    keyword = Stringmanip::to_upper(keyword);

    // read first word and check against keywordMap for action
    switch(keywordMap[keyword]){
      case NoMatch:
        throw PlatoExcept(keyword + " is not a valid keyword in config file");
        break;
      case keyAnalysis:
        // start new Analysis
        if(priorAnalysis){
          analyses.push_back(currAnalysis);
          optionVector.push_back(currListOptions);
          currAnalysis.clear();
        }
        priorAnalysis = true;
        break;
      case keyStart:
        {
          // check that option for loci management
          // has been set -- if not use default of passing list
          if(currListOptions.size() == currAnalysis.size())
            currListOptions.push_back(defaultOptions);

          FilterParams currFilter;
          ss >> currFilter.filter_name;
          currFilter.filter_name = Stringmanip::to_upper(currFilter.filter_name);
          read_params(currFilter, dataStream);
          currAnalysis.push_back(currFilter);
        }
        break;
      case keyEnd:
        throw PlatoExcept(keyword + " is unmatched by a START keyword");
        
        break;
      case keyDataset:
        {
          ss >> dataset_info;
        }
        break;
      case keyMissingValue:
        ss >> missingValue;
        if(missingValue >= 0 && missingValue <=2){
          throw PlatoExcept(keyword + " must be less than 0 or greater than 2");
        }
      case keyOut:
        {
          ss >> outputParam;
        }
      case keyRandSeed:
        { 
          ss >> randSeed;
        }
        break;
      case keyLocCombo:
        add_loci_options(ss, currListOptions); 
        break; 
      case keyDatasetType:
        ss >> dataset_type;
        break;
      case keyMapFile:
        ss >> loci_info_id;
        loci_info_available = true;
        break;
      case keyPseudoGenerate:
        ss >> pseudo_option;
        pseudo_option = Stringmanip::to_upper(pseudo_option);
        if(pseudo_option.compare("ON") == 0)
          pseudo_generate = true;
        break;
      case keyOutputName:
        ss >> output_name;
        break;
      default:
        throw PlatoExcept(keyword + " is not a valid keyword in config file");
    }; 
  }
  // should have one analysis left to push onto vector
  analyses.push_back(currAnalysis);
  optionVector.push_back(currListOptions);

  // establish output stream use input datafile with
  // .out appended to it

  dataStream.close();
}


///
/// Adds option for handling next filter in list 
/// @param ss stringstream with parameters 
/// @param optVect vector of ListOptions for handling loci
/// during transition from one filter to the next
/// @return
///
void ConfigFileReader::add_loci_options(stringstream & ss, 
  vector<ListOptions> & optVect){
  string keyword;
  ss >> keyword;
  keyword = Stringmanip::to_upper(keyword);
  ListOptions newOpt;

  // read first word and check against keywordMap for action
  switch(lociKeyMap[keyword]){
    case NoMatch:
      throw PlatoExcept(keyword + " is not a valid keyword in config file");
      break;
    case addLoci:
      newOpt.op = AddToList; 
      {
      unsigned int sizeCombo;
      do{
        ss >> sizeCombo; 
        newOpt.lociCombinations.push_back(sizeCombo);
      }while(!ss.eof());
      }
      if(newOpt.lociCombinations.size() > 1 && (newOpt.lociCombinations[newOpt.lociCombinations.size()-1] ==
        newOpt.lociCombinations[newOpt.lociCombinations.size()-2])){
          newOpt.lociCombinations.pop_back();
        }
      break;
    case removeLoci:
      newOpt.op = RemoveFromList;      
      {
      unsigned int sizeCombo;
      do{
        ss >> sizeCombo;
        newOpt.lociCombinations.push_back(sizeCombo);      
      }while(!ss.eof());
      }
      break;
    case passLoci:
      newOpt.op = PassList;
      break;
    default:
      throw PlatoExcept(keyword + " is not a valid keyword in config file");
  };

  optVect.push_back(newOpt);
}


///
/// Parses parameters for indicated filter and creates sub filters as needed
/// All subfilters must be nested within the filter using it
/// START LD
/// SUB CHI  --> start of chi square sub-filter
/// CHITYPE  SINGLE
/// THRESHOLD .05
/// END CHI   --> end of chi square sub-filter
/// THRESHOLD .05  --> this parameter goes with LD filter
/// @param  currFilter current filter
/// @param  configFile input stream
/// @return  none
///
void ConfigFileReader::read_params(FilterParams & currFilter, ifstream & configFile){

  string currLine, keyword;
  bool done = false;
  std::map<string, string> params;
  
  while(!configFile.eof() && !done){
    getline(configFile, currLine);
    
    if(skip_line(currLine)){
      continue;
    }
    
    stringstream ss(currLine+" ");
    ss >> keyword;
    keyword = Stringmanip::to_upper(keyword);

    switch(keywordMap[keyword]){
      case keySub:
        {
          FilterParams subFilter;
          // start a sub filter
          ss >> subFilter.filter_name;
          subFilter.filter_name = Stringmanip::to_upper(subFilter.filter_name);
          read_params(subFilter, configFile);
          currFilter.sub_filters.push_back(subFilter);
          
          break;
        }
      case keyEnd:
        return; // done reading parameters for this filter/subfilter
      default:
        {
          // add to parameter map for this filter
          string holder, joined;
          bool start = true;
          ss >> holder;
          do{
            if(start == false){
              joined += " ";
            }
            joined += holder;
            start = false;
            ss >> holder;
          }while(!ss.eof());
          currFilter.filter_params[keyword] = joined;
        }
    }
  }
}

///
/// Checks whether to skip blank lines or comment lines.
/// @param currLine string to check as comment or blank
/// @return true if skipping
///
bool ConfigFileReader::skip_line(string currLine){
    
     // skip blank lines
    if(currLine.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") == string::npos){
      return true;
    }

    int character = 0;
    // skip comments -- begin with '#'
    while(currLine[character] == ' ' || currLine[character] == '\t'){
      character++;
    }
    if(currLine[character] == '#'){
      return true;
    }
  return false;
}
}



