/*
 * InputProcess.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#include "InputProcess.h"

namespace po=boost::program_options;

namespace ProcessLib{

const std::string InputProcess::stepname = InputProcess::doRegister("load-data");

void InputProcess::process(Methods::DataSet& ds){
	setDataSet(ds);
	read();
}

po::options_description& InputProcess::appendOptions(po::options_description& opts){
	return DataLoader::addOptions(opts);
}

void InputProcess::parseOptions(const po::variables_map& vm){
	DataLoader::parseOptions(vm);
}

}
