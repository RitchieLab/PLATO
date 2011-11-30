//LogInfo.cpp

#include "LogInfo.h"

#include <fstream>
#include <iostream>

using namespace Methods;
namespace Filters{
///
/// Logs parameters in a combo generator.  Overwrites any
/// previous information.
/// @param generator ComboGenerator with parameters
/// @param filename name of file
///
void LogInfo::log_parameters(ComboGenerator& generator, string filename){

  fstream outfile;
  outfile.open(filename.c_str(), ios::out);
  if(!outfile.is_open()){
    throw PlatoExcept(filename + " unable to for logging combination parameters");
  }
  
  // output the parameters
  outfile << generator.param_j() << " " << generator.param_x() << " "
    << generator.param_kdec() << endl;
  // now need to output c parameter 
  // number of c values is ComboEnd + 3
  int total_c = generator.get_combo_max() + 3;
  outfile << total_c << endl;

  outfile << generator.param_c(0);

  for(int i = 1; i < total_c; i++)
    outfile << " " << generator.param_c(i);
  outfile << endl;
  
  outfile.close();
  
}
}

