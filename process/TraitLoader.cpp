#include "TraitLoader.h"

#include "data/DataSet.h"
#include "data/Sample.h"
#include "util/Logger.h"

namespace po = boost::program_options;

using std::string;
using std::stringstream;

using PLATO::Utility::Logger;
using PLATO::Data::Sample;
using PLATO::Data::DataSet;

namespace PLATO {
namespace ProcessLib {

const std::string TraitLoader::stepname = TraitLoader::doRegister("load-trait");

po::options_description& TraitLoader::appendOptions(
		po::options_description& opts) {

	po::options_description file_opts = FileLoader::getOptions();

	po::options_description trait_opts("Trait Loading Options");

	trait_opts.add_options()
		("ignore-error",po::bool_switch(&ignore_error),"Treat any conversion errors as missing")
		;

	opts.add(file_opts).add(trait_opts);

	return opts;
}

void TraitLoader::parseOptions(const po::variables_map& vm) {
	FileLoader::parseOptions(vm);
}

void TraitLoader::process(DataSet& ds) {
	load(ds);
}

void TraitLoader::processEntry(DataSet& ds, const string& name, Sample& s, const string& value){
	float val;

	if (!(stringstream(value) >> val)) {
		Logger::log_err("Unparseable value for '"
				+ name + "' for sample " + s.getID(),
				!ignore_error);
	} else {
		ds.addTrait(name, &s, val);
	}
}

}
}

