#ifndef PROCESSLIB_TRAITMISSING_FILTER_H   //change to actual process name
#define PROCESSLIB_TRAITMISSING_FILTER_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class TraitMissingFilter : public ProcessImpl<TraitMissingFilter>{
private:
	const static std::string stepname;

public:
	TraitMissingFilter() : ProcessImpl<TraitMissingFilter>(stepname, "Filter the traits based on missingness rate") {};
	virtual ~TraitMissingFilter(){};

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
