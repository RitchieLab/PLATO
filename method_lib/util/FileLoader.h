/*
 * FileLoader.h
 *
 *  Created on: Sep 22, 2014
 *      Author: jrw32
 */

#ifndef UTILITY_FILELOADER_H
#define UTILITY_FILELOADER_H

#include <string>
#include <vector>
#include <boost/program_options.hpp>

namespace PLATO{

namespace Data{
class DataSet;
class Sample;
}

namespace Utility{


class FileLoader {
public:
	FileLoader() {}
	virtual ~FileLoader() {}

protected:
	boost::program_options::options_description getOptions();
	void parseOptions(const boost::program_options::variables_map& vm);


	void load(Data::DataSet& ds);
	virtual void processEntry(Data::DataSet& ds, const std::string& name, Data::Sample& s, const std::string& value) = 0;

private:
	static void readSampleList(const std::vector<std::string>& in_list, std::set<std::string>& out_set);
	static void readSampleFile(const std::vector<std::string>& in_list, std::set<std::string>& out_set);
	static void addSampleToSet(const std::string& samp, std::set<std::string>& out_set);
	bool filterSample(const std::string& id, const std::string& fid="") const;

protected:
	std::vector<std::string> trait_fns;
	std::string missing_val;

	bool no_fid;
	bool extra_samples;
	bool dummy_samples;
	bool require_complete;

private:
	std::vector<std::string> incl_sample_str;
	std::vector<std::string> excl_sample_str;
	std::vector<std::string> incl_sample_fns;
	std::vector<std::string> excl_sample_fns;
	std::set<std::string> incl_sample_set;
	std::set<std::string> excl_sample_set;
	static const std::string sampl_field_sep;
};

}
}

#endif /* FILELOADER_H_ */
