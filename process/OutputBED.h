/*
 * OutputBED.h
 *
 *  Created on: Dec 3, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_OUTPUTBED_H
#define PROCESSLIB_OUTPUTBED_H

#include "Process.h"
#include "util/OutputPLINK.h"

namespace ProcessLib {

class OutputBED : public ProcessImpl<OutputBED>, private Methods::OutputPLINK {
private:
	const static std::string stepname;

public:
	OutputBED() : ProcessImpl<OutputBED>(stepname, "Output dataset in BED/BIM/FAM format"){}
	virtual ~OutputBED() {}

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:

	unsigned char getBinaryGeno(const Methods::Sample& s, const Methods::Marker& m) const;

private:
	std::string base_fn;
	std::string bed_fn;
	std::string bim_fn;
	std::string fam_fn;

	bool _ind_major;
};

}

#endif /* OUTPUTBED_H_ */
