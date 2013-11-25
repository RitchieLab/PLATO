#include "Process.h"

namespace po=boost::program_options;

Process::Process() : data_set(NULL), overwrite(false), order(0),
		orig_num_markers(0), _MARKERLIST_(false), _STRATIFY_(false) {
	options.setCovarMissing(Methods::opts::_COVAR_MISSING_);
	options.setTraitMissing(Methods::opts::_TRAIT_MISSING_);
}

po::options_description& Process::addOptions(po::options_description& opts){
	opt_ptr = &opts;
	return appendOptions(opts);
}

void Process::run(Methods::DataSet* ds){
	process(ds);
	PrintSummary();
}


