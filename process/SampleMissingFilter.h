#ifndef PROCESSLIB_SAMPLEMISSING_FILTER_H   //change to actual process name
#define PROCESSLIB_SAMPLEMISSING_FILTER_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class SampleMissingFilter : public ProcessImpl<SampleMissingFilter>{
private:
	const static std::string stepname;

public:
	SampleMissingFilter() : ProcessImpl<SampleMissingFilter>(stepname, "Filter the dataset based on sample call rate") {};
	virtual ~SampleMissingFilter(){};

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
