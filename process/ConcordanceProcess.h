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

namespace PLATO{
namespace ProcessLib {

class ConcordanceProcess : public ProcessImpl<ConcordanceProcess>, private Utility::DataLoader {

private:
	const static std::string stepname;

public:
	ConcordanceProcess() : ProcessImpl<ConcordanceProcess>(stepname, "Command to check concordance") {};
	virtual ~ConcordanceProcess(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

};

}
}

#endif /* ConcordanceProcess_H_ */
