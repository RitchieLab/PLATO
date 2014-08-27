/*
 * ConcordanceProcess.h
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_CONCORDNACEPROCESS_H
#define PROCESSLIB_CONCORDANCEPROCESS_H

#include "Process.h"
#include "util/DataLoader.h"

#include <iosfwd>

namespace PLATO{
namespace ProcessLib {

class ConcordanceProcess : public ProcessImpl<ConcordanceProcess>, public Utility::DataLoader {

private:
	const static std::string stepname;

public:
	ConcordanceProcess() : ProcessImpl<ConcordanceProcess>(stepname, "Command to check concordance") {};
	virtual ~ConcordanceProcess(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	// check to see if the genotypes match at a given point
	// return value & 1 gives if they did not match
	// return value & 2 gives if none were missing
	unsigned char match(std::ofstream* error_f, const Data::Marker* m, const Data::Sample* s,
			   const Data::Marker* a_m, const Data::Sample* a_s, bool phased) const;

	void outputMarkerMismatch(std::ofstream* f, const Data::Marker* m, bool in_base) const;
	void outputSampleMismatch(std::ofstream* f, const Data::Sample* s, bool in_base) const;

	void printMarker(std::ofstream& f, const Data::Marker& m) const;
	void printSample(std::ofstream& f, const Data::Sample& s, bool extra_data=true) const;

private:
	std::string _prefix;
	std::string _mmismatch_fn;
	std::string _smismatch_fn;
	std::string _error_fn;
	std::string _sample_fn;
	std::string _marker_fn;
	std::string _ext;
	std::string _sep;

	bool _incl_missing;
};

}
}

#endif /* ConcordanceProcess_H_ */
