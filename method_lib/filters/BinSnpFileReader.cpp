//BinSnpFileReader.cpp

#include "BinSnpFileReader.h"
#include "Stringmanip.h"

using namespace std;
namespace Filters{

///
/// Opens file appropriately for this SNP file reader
/// @param filename 
/// @return true when file opened successfully
///
bool BinSnpFileReader::open_file(string filename){
  filereader.open(filename.c_str(), ifstream::binary);
  // need to skip ahead past the header information in the file
  if(filereader.is_open()){
    filereader.seekg(2 * sizeof(unsigned int));
  }
  return filereader.is_open();
}


///
/// Returns set of models from SNP file
/// @param combos 2-D vector of strings with SNP names
/// @return True when end of file reached
///
bool BinSnpFileReader::get_combinations(vector<vector<string> >& combos, int nModels){
  // empty combinations
  combos.clear();
  
  bool endFile = false;
  vector<string> model(2,"");
  unsigned int rsone, rstwo;
  float implication;
  
  for(int combo=0; combo<nModels; combo++){
    if(filereader.eof()){
      endFile=true;
      break;
    }
  
    filereader.read( (char *)&rsone, sizeof(unsigned int));
    filereader.read( (char *)&rstwo, sizeof(unsigned int));
    filereader.read( (char *)&implication, sizeof(float));
    model[0] = "rs"+ Stringmanip::itos(rsone);
    model[1] = "rs"+ Stringmanip::itos(rstwo);
    combos.push_back(model);
  }
  
  return endFile;
  
}

}
