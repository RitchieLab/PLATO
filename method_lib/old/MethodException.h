//MethodException.h

#ifndef __METHODEXCEPT_H__
#define __METHODEXCEPT_H__

#include <exception>
#include <string>
using namespace std;

///
/// Class thrown for exception in data handling  <br>
/// Error messages are set by the creating class
///


/// Exception class
namespace Methods{
class MethodException: public std::exception{
        public:
          MethodException() throw(){};
          MethodException(string message){error=message;}
          ~MethodException() throw(){};
          virtual const char* what() const throw(){return error.c_str();}

        private:
          string error;
};
};
#endif

