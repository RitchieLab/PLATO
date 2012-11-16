//BioBinaryReader.cpp

#include "BioBinaryReader.h"
#include <sstream>
#include "Stringmanip.h"
#include "FilterExcept.h"
#include <iostream>
#include "BinSnpFileReader.h"
#include "TextSnpFileReader.h"

using namespace std;

namespace Filters{

///
/// Constructor
///
BioBinaryReader::BioBinaryReader(){
  snpreader = NULL;
}


BioBinaryReader::BioBinaryReader(string filename){
  snpreader = NULL;
  SetFilename(filename);
}


///
/// Destructor
///
BioBinaryReader::~BioBinaryReader(){
    snpreader->close_file();
    if(snpreader != NULL)
      delete snpreader;
}


///
/// Sets filename and opens file for reading
/// @param filename
///
void BioBinaryReader::SetFilename(string filename){

  // determine type of file to read and create appropriate
  // type of SNP File Reader to handler it
  
  // if it is binary the first unsigned int will equal 0
  // otherwise it will be some other number that forms part of
  // the count in the text file
  ifstream reader;
  unsigned int binaryCheck;
  reader.open(filename.c_str(), ios::binary);
  if(!reader.is_open()){
    throw FilterExcept("Unable to open bio filter file " + filename);
  }
  reader.read((char*)&binaryCheck, sizeof(unsigned int));
  reader.close();
  
  if(binaryCheck==0){
    snpreader = new BinSnpFileReader;
  }
  else{
    snpreader = new TextSnpFileReader;
  }
  
  snpreader->open_file(filename);
  
}


///
/// Gets next set of models from binary file and puts in linked list
/// @param combos
/// @return false when all combinations returned
///
bool BioBinaryReader::get_combinations(vector<vector<string> >& combos){
  
  return snpreader->get_combinations(combos, 30000);
  
}
}




