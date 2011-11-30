//BinSnpFileReader.h

#ifndef __BINSNPFILEREADER_H__
#define __BINSNPFILEREADER_H__

#include "SnpFileReader.h"

namespace Filters{
/// Abstract base class for objects that read different formats from biofilter application
class BinSnpFileReader: public SnpFileReader{

  public:

    bool open_file(std::string filename);

    bool get_combinations(std::vector<std::vector<std::string> >& combos, int nModels=30000);
    
};

}

#endif
