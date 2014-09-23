#ifndef PROCESSLIB_CATEGORICALLOADER_H
#define PROCESSLIB_CATEGORICALLOADER_H

#include "Process.h"
#include "util/FileLoader.h"

#include <map>
#include <string>

namespace PLATO{
namespace ProcessLib {

class CategoricalLoader : public ProcessImpl<CategoricalLoader>, public Utility::FileLoader {

private:
	const static std::string stepname;

public:
	CategoricalLoader() : ProcessImpl<CategoricalLoader>(stepname, "Command to load categorical covariate data into PLATO"),
		Utility::FileLoader(){};
	virtual ~CategoricalLoader(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

	virtual void processEntry(Data::DataSet& ds, const std::string& name, Data::Sample& s, const std::string& value);

private:
	std::map<std::string, std::map<std::string, unsigned char> > _name_val_map;
};

}
}

#endif /* INPUTPROCESS_H_ */
