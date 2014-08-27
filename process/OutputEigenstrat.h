#ifndef PROCESSLIB_OUTPUTEIGENSTRAT_H
#define PROCESSLIB_OUTPUTEIGENSTRAT_H

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class OutputEigenstrat : public ProcessImpl<OutputEigenstrat>{
private:
	const static std::string stepname;

	const static unsigned int missing_val = 9;

public:
	OutputEigenstrat() : ProcessImpl<OutputEigenstrat>(stepname, "Output in Eigenstrat format") {};
	virtual ~OutputEigenstrat(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void PrintSummary();
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string _prefix;
	std::string _geno_fn;
	std::string _snp_fn;
	std::string _indiv_fn;

};

}
}

#endif
