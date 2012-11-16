//FilterExcept.cpp

#include "FilterExcept.h"
namespace Filters{
const char* FilterExcept::what() const throw()
{
  return error.c_str();
}

FilterExcept::FilterExcept(string message){
  error = message;
}
}
