#ifndef PROCESSLIB_MARKERMISSING_FILTER_H   //change to actual process name
#define PROCESSLIB_MARKERMISSING_FILTER_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class MarkerMissingFilter : public ProcessImpl<MarkerMissingFilter>{
private:
	const static std::string stepname;

public:
	MarkerMissingFilter() : ProcessImpl<MarkerMissingFilter>(stepname, "Filter the dataset based on marker call rate") {};
	virtual ~MarkerMissingFilter(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void PrintSummary();
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:

	double _thresh;
	int _n_filtered;

};

}
}

#endif
