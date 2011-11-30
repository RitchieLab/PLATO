//BioBinaryReader.h

#ifndef __BIOBINARYREADER_H__
#define __BIOBINARYREADER_H__

#include <fstream>
#include <string>
#include <vector>
#include "SnpFileReader.h"

namespace Filters{
/// Get combinations from specified tables in database
class BioBinaryReader{

  public:
    /// Contructor
    BioBinaryReader();
    
    /// Alternative constructor
    BioBinaryReader(std::string filename);
    
    /// Destructor
    ~BioBinaryReader();
 
    /// Retrieve combinations -- return false when all combinations returned
    bool get_combinations(std::vector<std::vector<std::string> >& combos);
    
    void SetFilename(std::string filename);
        
  private:
    
    SnpFileReader* snpreader;

};

}

#endif
