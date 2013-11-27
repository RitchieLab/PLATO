/*
 * BatchProcess.h
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_BATCHPROCESS_H
#define PROCESSLIB_BATCHPROCESS_H

#include "Process.h"

#include <string>

namespace ProcessLib {

class BatchProcess : public ProcessImpl<BatchProcess> {
private:
	const static std::string stepname;

public:
	BatchProcess() : ProcessImpl<BatchProcess>(stepname, "Execute commands from a file") {}

	virtual void parseOptions(const boost::program_options::variables_map& vm) {}

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string batch_fn;

};

}

#endif /* BATCHPROCESS_H_ */
