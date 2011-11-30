//TextSnpFileReader.cpp

#include "TextSnpFileReader.h"

#include "Stringmanip.h"

using namespace std;

namespace Filters{

///
/// Opens file appropriately for this SNP file reader
/// @param filename 
/// @return true when file opened successfully
///
bool TextSnpFileReader::open_file(string filename){
  filereader.open(filename.c_str(), ifstream::in);
  // need to skip ahead past the header information in the file
  string totalModels;
  if(filereader.is_open()){
    filereader >> totalModels;
  }
  return filereader.is_open();
}


///
/// Returns set of models from SNP file
/// @param combos 2-D vector of strings with SNP names
/// @return True when end of file reached
///
bool TextSnpFileReader::get_combinations(vector<vector<string> >& combos, int nModels){
  // empty combinations
  combos.clear();
  
  bool endFile = false;
  vector<string> model(2,"");

  string rsone, rstwo, implication;
  
  for(int combo=0; combo<nModels; combo++){
    if(filereader.eof()){
      endFile=true;
      break;
    }
  
    filereader >> rsone >> rstwo >> implication;
    model[0] = "rs" + rsone;
    model[1] = "rs" + rstwo;
    
    combos.push_back(model);
  }
  
  return endFile;
  
}

}

