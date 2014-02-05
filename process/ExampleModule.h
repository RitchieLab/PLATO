#ifndef PROCESSLIB_EXAMPLEMODULE_H   //change to actual process name
#define PROCESSLIB_EXAMPLEMODULE_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class ExampleModule : public ProcessImpl<ExampleModule>{
private:
	const static std::string stepname;

public:
	ExampleModule() : ProcessImpl<ExampleModule>(stepname, "Example Command - please remove in released version") {};
	virtual ~ExampleModule(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void PrintSummary();
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string arg_string;
	bool arg_bool;

};

}
}

#endif
