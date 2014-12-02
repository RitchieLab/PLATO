/*
 * FileLoader.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: jrw32
 */

#include "FileLoader.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "util/Logger.h"
#include "data/DataSet.h"
#include "data/Sample.h"

using std::string;
using std::vector;
using std::stringstream;
using std::ifstream;
using std::set;

using boost::algorithm::split;
using boost::algorithm::trim;
using boost::algorithm::is_any_of;

using PLATO::Utility::Logger;
using PLATO::Data::DataSet;
using PLATO::Data::Sample;

namespace po=boost::program_options;

namespace PLATO{
namespace Utility{

po::options_description FileLoader::getOptions(){

	po::options_description file_opts("File Loading Options");

	file_opts.add_options()("file", po::value<vector<string> >(&trait_fns)->composing(),"Trait file to load")
			("missing", po::value<string>(&missing_val),"Missing value")
			("no-fid", po::bool_switch(&no_fid),"Trait file has no FamID column")
			("extra-samples", po::bool_switch(&extra_samples),"Ignore any samples that cannot be mapped to existing data")
			("dummy-samples", po::bool_switch(&dummy_samples),"Create samples for any that cannot be mapped to existing data")
			("require-complete", po::bool_switch(&require_complete), "Require trait data for all (enabled) samples loaded so far");

	return file_opts;
}

void FileLoader::load(DataSet& ds){
	vector<string>::const_iterator fn_itr = trait_fns.begin();

	while (fn_itr != trait_fns.end()) {

		ifstream input((*fn_itr).c_str());

		if (!input.is_open()) {
			Logger::log_err("Error opening trait file: " + (*fn_itr), true);
		}

		string line;
		getline(input, line);
		trim(line);
		vector<string> headers;
		split(headers, line, is_any_of(" \n\t"), boost::token_compress_on);

		vector<string> values;
		values.reserve(headers.size());
		int lineno = 1;
		set<const Sample*> enabled_samples(ds.beginSample(), ds.endSample());
		unsigned int n_extra = 0;

		while (getline(input, line)) {
			++lineno;
			trim(line);
			split(values, line, is_any_of(" \n\t"), boost::token_compress_on);

			Sample* s = 0;
			if (no_fid) {
				s = ds.getSample(values[0]);
			} else {
				s = ds.getSample(values[0], values[1]);
			}

			if (!s) {
				n_extra++;
				if(!(extra_samples || dummy_samples)){
					Logger::log_err(string("Extra sample found on line ") + boost::lexical_cast<string>(lineno), true);
				}
				if (dummy_samples) {
					if (no_fid) {
						s = ds.addSample(values[0]);
					} else {
						s = ds.addSample(values[0], values[1]);
					}
				}

			} else {
				enabled_samples.erase(s);
			}

			if (s) {
				for (unsigned int i = (1 + (!no_fid)); i < std::min(
						values.size(), headers.size()); i++) {
					if (values[i] != missing_val) {
						processEntry(ds, headers[i], *s, values[i]);
					}
				}
			}
		}

		input.close();

		if(n_extra > 0){
			Logger::log_err("INFO: " + boost::lexical_cast<string>(n_extra) + " extra samples found in " + *fn_itr);
		}

		if(enabled_samples.size() > 0){
			Logger::log_err("Trait data not given in " + *fn_itr + " for some previously loaded samples.", require_complete);
		}
		++fn_itr;
	}
}

}
}
