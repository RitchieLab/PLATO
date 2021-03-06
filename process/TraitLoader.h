#ifndef PROCESSLIB_TRAITLOADER_H
#define PROCESSLIB_TRAITLOADER_H

#include "Process.h"
#include "util/FileLoader.h"

#include <string>

namespace PLATO{
namespace ProcessLib {

class TraitLoader : public ProcessImpl<TraitLoader>, public Utility::FileLoader {

private:
	const static std::string stepname;

public:
	TraitLoader() : ProcessImpl<TraitLoader>(stepname, "Command to load numeric trait/covariate data into PLATO"), Utility::FileLoader() {};
	virtual ~TraitLoader(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

	virtual void processEntry(Data::DataSet& ds, const std::string& name, Data::Sample& s, const std::string& value);

private:
	bool ignore_error;

};

}
}

#endif /* INPUTPROCESS_H_ */
