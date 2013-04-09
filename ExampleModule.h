#ifndef EXAMPLEMODULE_H   //change to actual process name
#define EXAMPLEMODULE_H   //change to actual process name

#include <string>

#include "Process.h"

class ExampleModule : public ProcessImpl<ExampleModule>{
private:
	const static std::string stepname;

public:
	ExampleModule(){};
	virtual ~ExampleModule(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);
};

#endif
