//OutputWriterFactory.h
#ifndef __OUTPUTWRITERFACTORY_H__
#define __OUTPUTWRITERFACTORY_H__

#include <map>
#include <string.h>
#include <string>

using namespace std;

///
/// Creates OutputWriter for use in output of results
///
namespace Filters{

class OutputWriter;

/// Creates OutputWriter
class OutputWriterFactory{
  public:
    /// Creates filter and returns.  Throws exception if 
    /// no filter matches that name
    static OutputWriter * create_writer(std::string outputName);
    
    /// Enumeration for filter types with a new type added for each new filter in system
    enum OutputWriterType{
      NoMatch,
      /// Enum for Database storage of results
      DBOutputType, 
      /// Enum for Text file storage of results
      TextFileOutputType,
      /// Enum for Will's database testing
      LRDBOutputType
    }; 

  private: 
    /// Sets map for use in Filter creation
    static void SetFilterMap();
    
    static std::map<std::string, OutputWriterType> OutputMap;

};
}

#endif

