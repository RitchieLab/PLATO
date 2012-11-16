//FilterExcept.h

#ifndef __FILTEREXCEPT_H__
#define __FILTEREXCEPT_H__

#include <exception>
#include <string>
using namespace std;

///
/// Class thrown for exception in data handling  <br>
/// Error messages are set by the creating class
///

namespace Filters{

/// Exception class
class FilterExcept: public std::exception{
        public:
          FilterExcept() throw();
          FilterExcept(string message);
          ~FilterExcept() throw(){};
          virtual const char* what() const throw();

        private:
          string error;
};
}

#endif

