#ifndef PROCESSLIB_TRAITLOADER_H
#define PROCESSLIB_TRAITLOADER_H

#include "Process.h"

#include <string>
#include <vector>
namespace PLATO{
namespace ProcessLib {

class TraitLoader : public ProcessImpl<TraitLoader> {

private:
	const static std::string stepname;

public:
	TraitLoader() : ProcessImpl<TraitLoader>(stepname, "Command to load numeric trait/covariate data into PLATO") {};
	virtual ~TraitLoader(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::vector<std::string> trait_fns;
	std::string missing_val;

	bool no_fid;
	bool ignore_error;
	bool extra_samples;
};

}
}

#endif /* INPUTPROCESS_H_ */
