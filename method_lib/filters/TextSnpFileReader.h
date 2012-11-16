//TextSnpFileReader.h

#ifndef __TEXTSNPFILEREADER_H__
#define __TEXTSNPFILEREADER_H__

#include "SnpFileReader.h"

namespace Filters{
/// Abstract base class for objects that read different formats from biofilter application
class TextSnpFileReader: public SnpFileReader{

  public:

    bool open_file(std::string filename);

    bool get_combinations(std::vector<std::vector<std::string> >& combos, int nModels=30000);
    
};

}

#endif
