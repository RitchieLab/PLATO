//LogInfo.h

#ifndef __LOGINFO_H__
#define __LOGINFO_H__

#include "ComboGenerator.h"
#include "PlatoExcept.h"

using namespace std;

///
/// Logs run information for later use
///

namespace Filters{

/// Contains filters for running analysis
class LogInfo{
  
  public:
    /// logs parameters in a combo generator
    void log_parameters(Methods::ComboGenerator& generator, string filename);
  
  private:
    
};

}


#endif
