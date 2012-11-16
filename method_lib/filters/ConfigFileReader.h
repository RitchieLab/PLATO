// ConfigFileReader.h

#ifndef __CONFIGFILEREADER_H__
#define __CONFIGFILEREADER_H__

#include "ConfigReader.h"
#include "Stringmanip.h"
#include <fstream>
#include <iostream>

///
/// Reads configuration file for plato
/// and stores parameters for access by
/// other objects
///
namespace Filters{
/// Reads configuration file and stores parameters
class ConfigFileReader: public ConfigReader{

  public:
  
    /// Destructor
    virtual ~ConfigFileReader(){}
  
    /// Default constructor
    ConfigFileReader(){initialize_keywords();}
  
    /// Alternate constructor
    ConfigFileReader(string config_file);
  
    /// Read config file and set filter parameters
    virtual void read_config(string configfile);   
   
  private:
    
    /// Reads parameters for a single filter
    void read_params(FilterParams & currFilter, ifstream & configFile);
    
    /// Skips line in configuration file
    bool skip_line(string currLine);

    /// Initializes keywords for 
    void initialize_keywords();

    /// Adds loci handling option to the vector
    void add_loci_options(stringstream & ss, vector<ListOptions> & optVect);

    enum configKeywords{
      keyNoMatch, 
      keyAnalysis,
      keyStart,
      keyEnd,
      keySub,
      keyDataset,
      keyOut,
      keyRandSeed,
      keyMissingValue,
      keyLocCombo,
      keyDatasetType,
      keyMapFile,
      keyPseudoGenerate,
      keyOutputName
    };
    
    enum lociKeywords{
      NoMatch,
      addLoci,
      removeLoci,
      passLoci
    };
   
    std::map<string, configKeywords> keywordMap; 
    std::map<string, lociKeywords> lociKeyMap;
};
}
#endif
