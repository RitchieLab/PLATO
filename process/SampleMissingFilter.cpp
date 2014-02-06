#include "SampleMissingFilter.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"

#include "util/Logger.h"
#include "util/Missing.h"

#include <boost/lexical_cast.hpp>

using std::string;

using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string SampleMissingFilter::stepname = SampleMissingFilter::doRegister("sample-call-filter");

//PrintSummary()
//used to output results after processing the data
void SampleMissingFilter::PrintSummary(){
	//output for process goes here
	//
	Utility::Logger::log("Sample Call Filter: Dropped " + boost::lexical_cast<string>(_n_filtered) + " Samples.\n");
}

//process()
//main method to get the process going and to the work
void SampleMissingFilter::process(DataSet& ds){

	_n_filtered = 0;
	DataSet::sample_iterator si = ds.beginSample();
	while(si != ds.endSample()){
		if(1 - Utility::Missing::sampleMissing(ds, **si) < _thresh){
			(*si)->setEnabled(false);
		}

		++si;
	}
}

po::options_description& SampleMissingFilter::appendOptions(po::options_description& opts){
	po::options_description subopts("Sample Call Rate Options");

	subopts.add_options()
		("threshold", po::value<double>(&_thresh),"Drop all samples with a call rate below this threshold")
		;

	opts.add(subopts);
	return opts;
}

void SampleMissingFilter::parseOptions(const po::variables_map& vm){
	// no need to do anything here!
}

}
}
