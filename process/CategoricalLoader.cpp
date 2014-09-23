#include "CategoricalLoader.h"


#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "util/Logger.h"

#include "data/DataSet.h"
#include "data/Sample.h"

namespace po = boost::program_options;

using std::string;
using std::map;

using PLATO::Utility::Logger;
using PLATO::Data::Sample;
using PLATO::Data::DataSet;

namespace PLATO {
namespace ProcessLib {

const std::string CategoricalLoader::stepname = CategoricalLoader::doRegister("load-categorical");

po::options_description& CategoricalLoader::appendOptions(
		po::options_description& opts) {

	po::options_description file_opts = FileLoader::getOptions();
	po::options_description trait_opts("Categorical Loading Options");


	opts.add(file_opts).add(trait_opts);

	return opts;
}

void CategoricalLoader::parseOptions(const po::variables_map& vm) {
	FileLoader::parseOptions(vm);
}

void CategoricalLoader::process(DataSet& ds) {
	load(ds);
}

void CategoricalLoader::processEntry(DataSet& ds, const string& name, Sample& s, const string& value){
	map<string, unsigned char>& val_map = _name_val_map[name];
	map<string, unsigned char>::const_iterator vitr = val_map.find(value);
	unsigned char curr_val = val_map.size();

	if(vitr == val_map.end()){
		val_map[value] = curr_val;
	} else {
		curr_val = (*vitr).second;
	}

	ds.addCategorical(name, &s, curr_val);
}



}
}

