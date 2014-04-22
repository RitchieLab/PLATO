#include "Process.h"

#include <boost/lexical_cast.hpp>
#include <boost/date_time.hpp>

#include "util/Logger.h"

namespace po=boost::program_options;

using PLATO::Utility::Logger;
using std::string;

namespace PLATO{

po::options_description& Process::addOptions(po::options_description& opts){
	opt_ptr = &opts;
	return appendOptions(opts);
}

void Process::run(PLATO::Data::DataSet& ds){
	Logger::log(boost::lexical_cast<string>(boost::posix_time::second_clock::local_time())
			+ ": Starting " + getName());
	process(ds);
	PrintSummary();
	Logger::log(boost::lexical_cast<string>(boost::posix_time::second_clock::local_time())
			+ ": Finished " + getName());
}

}


