#ifndef PROCESSLIB_ALLELEFREQ_FILTER_H   //change to actual process name
#define PROCESSLIB_ALLELEFREQ_FILTER_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class AlleleFreqFilter : public ProcessImpl<AlleleFreqFilter>{
private:
	const static std::string stepname;

public:
	AlleleFreqFilter() : ProcessImpl<AlleleFreqFilter>(stepname, "Filter the dataset based on minor allele frequency") {};
	virtual ~AlleleFreqFilter(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void PrintSummary();
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:

	double _min_thresh;
	double _max_thresh;
	int _n_filtered;
	bool _noconvert;

};

}
}

#endif
