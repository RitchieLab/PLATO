/*
 * InputProcess.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#include "InputProcess.h"

using PLATO::Utility::DataLoader;
using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string InputProcess::stepname = InputProcess::doRegister("load-data");

void InputProcess::process(DataSet& ds){
	read(ds);
}

po::options_description& InputProcess::appendOptions(po::options_description& opts){
	return DataLoader::addOptions(opts);
}

void InputProcess::parseOptions(const po::variables_map& vm){
	DataLoader::parseOptions(vm);
}

}
}
