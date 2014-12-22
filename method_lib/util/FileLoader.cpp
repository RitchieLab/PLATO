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

// Define the "Group Separator" string as the separation between FID and IID!
// Do this b/c it's legal ASCII, but darn near impossible to include on the cmd line
const string FileLoader::sampl_field_sep = "\x1d";

po::options_description FileLoader::getOptions(){

	po::options_description file_opts("File Loading Options");

	file_opts.add_options()("file", po::value<vector<string> >(&trait_fns)->composing(),"Trait file to load")
			("missing", po::value<string>(&missing_val),"Missing value")
			("no-fid", po::bool_switch(&no_fid),"Trait file has no FamID column")
			("extra-samples", po::bool_switch(&extra_samples),"Ignore any samples that cannot be mapped to existing data")
			("dummy-samples", po::bool_switch(&dummy_samples),"Create samples for any that cannot be mapped to existing data")
			("require-complete", po::bool_switch(&require_complete), "Require trait data for all (enabled) samples loaded so far");

	po::options_description filter_opts("Pre-filtering Options");

	filter_opts.add_options()
			("incl-sample", po::value<vector<string> >(&incl_sample_str)->composing(), "Sample(s) to include")
			("excl-sample", po::value<vector<string> >(&excl_sample_str)->composing(), "Sample(s) to exclude")
			("incl-sample-fn", po::value<vector<string> >(&incl_sample_fns)->composing(), "File of sample(s) to include")
			("excl-sample-fn", po::value<vector<string> >(&excl_sample_fns)->composing(), "File of sample(s) to exclude")
			;

	return file_opts.add(filter_opts);
}

void FileLoader::parseOptions(const po::variables_map& vm) {

	if(!dummy_samples && !(incl_sample_str.empty() && excl_sample_str.empty())){
		Logger::log_err("WARNING: include/exclude sample list given without specifying --dummy-samples.  "
				"Ignoring the include/exclude list(s)");
	} else if(dummy_samples) {

		readSampleList(incl_sample_str, incl_sample_set);
		readSampleList(excl_sample_str, excl_sample_set);
		readSampleFile(incl_sample_fns, incl_sample_set);
		readSampleFile(excl_sample_fns, excl_sample_set);
	}
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
		unsigned int n_added = 0;

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
				++n_extra;
				if(!(extra_samples || dummy_samples)){
					Logger::log_err(string("Extra sample found on line ") + boost::lexical_cast<string>(lineno), true);
				}
				if (dummy_samples) {
					if (no_fid) {
						if(filterSample(values[0])){
							s = ds.addSample(values[0]);
							++n_added;
						}
					} else {
						if(filterSample(values[0], values[1])){
							s = ds.addSample(values[0], values[1]);
							++n_added;
						}
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
			Logger::log_err("INFO: " + boost::lexical_cast<string>(n_extra) +
					" extra samples found in " + *fn_itr + ", " +
					boost::lexical_cast<string>(n_added) + " samples created.");
		}

		if(enabled_samples.size() > 0){
			Logger::log_err("Trait data not given in " + *fn_itr + " for some previously loaded samples.", require_complete);
		}
		++fn_itr;
	}
}

bool FileLoader::filterSample(const string& id, const string& fid) const {
	string key = id + sampl_field_sep + (fid.empty() ? id : fid);
	return (incl_sample_set.empty() || incl_sample_set.count(key) > 0) &&
			(excl_sample_set.empty() || excl_sample_set.count(key) == 0);
}

void FileLoader::readSampleList(const vector<string>& in_list, set<string>& out_set){
	vector<string> all_samps;
	InputManager::parseInput(in_list, all_samps);

	for(unsigned int i=0; i<all_samps.size(); i++){
		addSampleToSet(all_samps[i], out_set);
	}
}

void FileLoader::readSampleFile(const vector<string>& fn_list, set<string>& out_set){
	string samp;
	for(unsigned int i=0; i<fn_list.size(); i++){
		ifstream f(fn_list[i].c_str());
		while(getline(f, samp)){
			boost::algorithm::trim(samp);
			addSampleToSet(samp, out_set);
		}
	}
}

void FileLoader::addSampleToSet(const string& samp, set<string>& out_set){
	stringstream ss(samp);
	string id, fid;
	if(ss >> id){
		fid = id;
		if(ss >> fid){
			Logger::log_err("WARNING: Sample '" + samp +
					"' has more than an 'FID IID', "
					"ignoring everything after 2nd space!");
		}
	}
	out_set.insert(id + sampl_field_sep + fid);
}

}
}
