//PlatoExcept.h
#ifndef __PLATOEXCEPT_H__
#define __PLATOEXCEPT_H__

#include <string>
#include <exception>
using namespace std;

///
/// Class thrown for exception in PLATO system <br>
/// Error messages are set by the creating class
///
namespace Filters{
/// Exception for plato
class PlatoExcept: public exception{
        public:
          PlatoExcept() throw();
          PlatoExcept(string message){error=message;}
          ~PlatoExcept()throw(){};
          virtual const char* what() const throw(){return error.c_str();}

        private:
          string error;
};
}

#endif

