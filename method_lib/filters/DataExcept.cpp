//DataExcept.cpp

#include "DataExcept.h"

namespace Filters{
const char* DataExcept::what() const throw()
{
  return error.c_str();
}

DataExcept::DataExcept(string message){
  error = message;
}
}


