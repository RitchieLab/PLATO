#include "TraitMissingFilter.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"

#include "util/Logger.h"
#include "util/Missing.h"

#include <boost/lexical_cast.hpp>

using std::string;

using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string TraitMissingFilter::stepname = TraitMissingFilter::doRegister("filter-trait-missing");

//PrintSummary()
//used to output results after processing the data
void TraitMissingFilter::PrintSummary(){
	//output for process goes here
	//
	Utility::Logger::log("Trait Missing Filter: Dropped " + boost::lexical_cast<string>(_n_filtered) + " Traits.\n");
}

//process()
//main method to get the process going and to the work
void TraitMissingFilter::process(DataSet& ds){

	_n_filtered = 0;
	DataSet::trait_iterator ti = ds.beginTrait();
	while(ti != ds.endTrait()){
		if(1 - Utility::Missing::traitMissing(ds, *ti) < _thresh){
			ds.setTraitEnabled(*ti, false);
		}

		++ti;
	}
}

po::options_description& TraitMissingFilter::appendOptions(po::options_description& opts){
	po::options_description subopts("Sample Call Rate Options");

	subopts.add_options()
		("threshold", po::value<double>(&_thresh)->default_value(0.99, "0.99"),"Drop all traits with a missing rate below this threshold")
		;

	opts.add(subopts);
	return opts;
}

void TraitMissingFilter::parseOptions(const po::variables_map& vm){
	// no need to do anything here!
}

}
}
