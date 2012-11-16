//SnpFileReader.h

#ifndef __SNPFILEREADER_H__
#define __SNPFILEREADER_H__

#include <fstream>
#include <string>
#include <vector>
#include "SnpFileReader.h"

namespace Filters{
/// Abstract base class for objects that read different formats from biofilter application
class SnpFileReader{

  public:

    virtual ~SnpFileReader(){}
    
    virtual bool get_combinations(std::vector<std::vector<std::string> >& combos, 
      int nModels=30000)=0;

    virtual bool open_file(std::string filename)=0;

    void close_file(){if(filereader.is_open()) filereader.close();}

  protected:
    std::ifstream filereader;
};

}

#endif
