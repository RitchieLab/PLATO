//ConfigReader.h

#ifndef __CONFIGREADER_H__
#define __CONFIGREADER_H__

#include "PlatoExcept.h"
#include "PlatoDefinitions.h"
#include "FilterDefinitions.h"

///
/// Reads configuration file for plato
/// and stores parameters for access by
/// other objects
///
namespace Filters{

/// Reads configuration file and stores parameters
class ConfigReader{

  public:
  
    virtual ~ConfigReader(){}
  
    /// Read config file and set filter parameters
    virtual void read_config(string configID)=0;
    
    /// Returns number of analysis sets in the configuration file
    unsigned int num_sets(){return analyses.size();}
    
    /// Returns parameters for the indicated set
    vector<FilterParams> get_params(int index){return analyses[index];}

    /// Returns options for the indicated set
    vector<ListOptions> get_loci_options(int index){return optionVector[index];}
    
    /// Returns dataset name
    string get_dataset_info(){return dataset_info;}
    
    /// Returns random seed
    unsigned int get_seed(){return randSeed;}
    
    /// Returns output parameter
    string get_output(){return outputParam;}

    /// Returns missing value when parsing datafile
    int get_missing_value(){return missingValue;}

    /// Returns dataset type
    string get_dataset_type(){return dataset_type;}
    
    /// Returns true if loci information is available
    bool loci_info(){return loci_info_available;}
    
    /// Returns loci information identifier
    string get_loci_info_id(){return loci_info_id;}
    
    /// Returns true if pseudocontrols should be generated for the dataset
    bool generate_pseudo_controls(){return pseudo_generate;}

    /// Returns true if output name was specified in the configuration file
    bool output_name_set(){return output_specified;}
    
    /// Returns output name
    string get_output_name(){return output_name;}

  protected:

    vector<vector<FilterParams> > analyses;
    unsigned int randSeed;
    int missingValue;
    string outputParam, dataset_info, dataset_type, loci_info_id, output_name;
    bool loci_info_available, pseudo_generate, output_specified;
    vector<vector<ListOptions> > optionVector;    

    
  private:

};
}
#endif

