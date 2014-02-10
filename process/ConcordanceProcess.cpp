/*
 * ConcordanceProcess.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#include "ConcordanceProcess.h"

#include "data/DataSet.h"

using PLATO::Utility::DataLoader;
using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string ConcordanceProcess::stepname = ConcordanceProcess::doRegister("concordance");

void ConcordanceProcess::process(DataSet& ds){
	DataSet alt_ds;
	read(alt_ds);

	// OK, now I can check the concordance between ds and alt_ds!

}

po::options_description& ConcordanceProcess::appendOptions(po::options_description& opts){
	return DataLoader::addOptions(opts);
}

void ConcordanceProcess::parseOptions(const po::variables_map& vm){
	DataLoader::parseOptions(vm);
}

}
}
