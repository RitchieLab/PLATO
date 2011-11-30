//DataExcept.h

#ifndef datahandling__DATAEXCEPT_H__
#define datahandling__DATAEXCEPT_H__

#include <string>
#include <exception>
using namespace std;


///
/// Class thrown for exception in data handling  <br>
/// Error messages are set by the creating class
///

namespace Filters{

/// Exception class
class DataExcept: public exception{
  public:
    DataExcept() throw();
    DataExcept(string message);
    ~DataExcept() throw(){};
    virtual const char* what() const throw();

  private:
    string error;

};

}

#endif
