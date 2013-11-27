/*
 * OutputPED.h
 *
 *  Created on: Nov 27, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_OUTPUTPED_H
#define PROCESSLIB_OUTPUTPED_H

#include "Process.h"

namespace ProcessLib {

class OutputPED : public ProcessImpl<OutputPED> {
private:
	const static std::string stepname;

public:
	OutputPED() : ProcessImpl<OutputPED>(stepname, "Output dataset in PED/MAP format"){}
	virtual ~OutputPED() {}

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string base_fn;
	std::string ped_fn;
	std::string map_fn;
};

}

#endif /* OUTPUTPED_H_ */
