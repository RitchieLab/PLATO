#include "Process.h"

namespace po=boost::program_options;

namespace PLATO{

po::options_description& Process::addOptions(po::options_description& opts){
	opt_ptr = &opts;
	return appendOptions(opts);
}

void Process::run(PLATO::Data::DataSet& ds){
	process(ds);
	PrintSummary();
}

}


