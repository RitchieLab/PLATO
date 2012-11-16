//FilterFactory.cpp
#include "OutputWriterFactory.h"
// #include "DBResults.h"
#include "FileResults.h"
// #include "BioDBResults.h"
namespace Filters{

std::map<std::string, OutputWriterFactory::OutputWriterType> OutputWriterFactory::OutputMap;

///
/// Function that creates a filter based on the name 
/// Used by system code in managing analysis pipeline
/// @param filterName name of filter
/// @return pointer to new Filter object
/// @throws PlatoExcept if not a valid filter name
///
OutputWriter * OutputWriterFactory::create_writer(std::string writerName){
  if(OutputMap.empty()){
    SetFilterMap();
  }
   
  OutputWriter * newWriter;
  switch(OutputMap[writerName]){
//     case DBOutputType:
//       newWriter = new DBResults;
//       break;
    case TextFileOutputType:
      newWriter = new FileResults;
      break;
//     case LRDBOutputType:
//       newWriter = new BioDBResults;
//       break;
    case NoMatch:
      throw PlatoExcept(writerName + " is not a valid output specifier.");
      break;
    default:
      throw PlatoExcept(writerName + " is not a valid output specifier.");
  };
  return newWriter;
}


///
/// Establishes the map for use in creating filters
/// @return 
///
void OutputWriterFactory::SetFilterMap(){
  OutputMap["DB"]=DBOutputType;
  OutputMap["TEXT"]=TextFileOutputType;
  OutputMap["LRWILL"]=LRDBOutputType;
}
}
