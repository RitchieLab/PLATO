/*
 * InputProcess.h
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#ifndef PROCESSLIB_INPUTPROCESS_H
#define PROCESSLIB_INPUTPROCESS_H

#include "Process.h"
#include "util/DataLoader.h"

namespace PLATO{
namespace ProcessLib {

class InputProcess : public ProcessImpl<InputProcess>, private Utility::DataLoader {

private:
	const static std::string stepname;

public:
	InputProcess() : ProcessImpl<InputProcess>(stepname, "Command to load data into PLATO"), DataLoader() {};
	virtual ~InputProcess(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

};

}
}

#endif /* INPUTPROCESS_H_ */
