#include "TraitLoader.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "util/Logger.h"

#include "data/DataSet.h"
#include "data/Sample.h"

namespace po=boost::program_options;

using std::string;
using std::vector;
using std::stringstream;
using std::ifstream;

using boost::algorithm::split;
using boost::algorithm::is_any_of;

using PLATO::Utility::Logger;
using PLATO::Data::Sample;
using PLATO::Data::DataSet;

namespace PLATO{
namespace ProcessLib {

const std::string TraitLoader::stepname =
		TraitLoader::doRegister("load-trait");

void TraitLoader::process(DataSet& ds) {
	vector<string>::const_iterator fn_itr = trait_fns.begin();

	while (fn_itr != trait_fns.end()) {

		ifstream input((*fn_itr).c_str());

		if (!input.is_open()) {
			throw std::invalid_argument("Error opening trait file: " + (*fn_itr));
		}

		string line;
		getline(input, line);
		vector<string> headers;
		split(headers, line, is_any_of(" \n\t"), boost::token_compress_on);

		vector<string> values;
		values.reserve(headers.size());
		int lineno = 1;
		while (getline(input, line)) {
			++lineno;
			split(values, line, is_any_of(" \n\t"), boost::token_compress_on);


			Sample* s = 0;
			if(no_fid){
				s = ds.getSample(values[0]);
			} else {
				s = ds.getSample(values[0], values[1]);
			}

			if(!s){
				string err("Extra sample found on line " +  static_cast<stringstream*>( &(stringstream() << lineno) )->str());
				Logger::log_err(err, !extra_samples);
			}else{
				for(unsigned int i=(1+(!no_fid)); i < std::min(values.size(), headers.size()); i++){
					float val;
					if(values[i] != missing_val && !(stringstream(values[i]) >> val)){
						string err("Unparseable value for '" + headers[i] +
							"' on line " +  static_cast<stringstream*>( &(stringstream() << lineno) )->str());
						Logger::log_err(err, !ignore_error);
					} else {
						ds.addTrait(headers[i], s, val);
					}
				}
			}
		}

		input.close();
		++fn_itr;
	}

}

po::options_description& TraitLoader::appendOptions(po::options_description& opts){
	po::options_description trait_opts("Trait Loading Options");

	trait_opts.add_options()
		("file", po::value<vector<string> >(&trait_fns)->composing(), "Trait file to load")
		("missing", po::value<string>(&missing_val), "Missing value")
		("no-fid", po::bool_switch(&no_fid), "If given, trait file has no FamID column")
		("ignore-error", po::bool_switch(&ignore_error), "If given, treat any conversion errors as missing")
		("extra-samples", po::bool_switch(&extra_samples), "If given, ignore any samples that cannot be mapped to existing data")
		;

	opts.add(trait_opts);

	return opts;
}

void TraitLoader::parseOptions(const po::variables_map& vm){
}

}
}

