/*
 * OutputPED.h
 *
 *  Created on: Nov 27, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_OUTPUTPED_H
#define PROCESSLIB_OUTPUTPED_H

#include "Process.h"
#include "util/OutputPLINK.h"

namespace ProcessLib {

class OutputTPED : public ProcessImpl<OutputTPED>, private Methods::OutputPLINK {
private:
	const static std::string stepname;

public:
	OutputTPED() : ProcessImpl<OutputTPED>(stepname, "Output dataset in transposed PED/MAP format"){}
	virtual ~OutputTPED() {}

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string base_fn;
	std::string tped_fn;
	std::string tfam_fn;
};

}

#endif /* OUTPUTPED_H_ */
