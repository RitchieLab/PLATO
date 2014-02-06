#include "MarkerMissingFilter.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"
#include "data/Marker.h"

#include "util/Logger.h"
#include "util/Missing.h"

#include <boost/lexical_cast.hpp>

using std::string;

using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string MarkerMissingFilter::stepname = MarkerMissingFilter::doRegister("marker-call-filter");

//PrintSummary()
//used to output results after processing the data
void MarkerMissingFilter::PrintSummary(){
	//output for process goes here
	//
	Utility::Logger::log("Marker Call Filter: Dropped " + boost::lexical_cast<string>(_n_filtered) + " Markers.\n");
}

//process()
//main method to get the process going and to the work
void MarkerMissingFilter::process(DataSet& ds){

	_n_filtered = 0;
	DataSet::marker_iterator mi = ds.beginMarker();
	while(mi != ds.endMarker()){
		if(1 - Utility::Missing::markerMissing(ds, **mi) < _thresh){
			(*mi)->setEnabled(false);
		}

		++mi;
	}
}

po::options_description& MarkerMissingFilter::appendOptions(po::options_description& opts){
	po::options_description subopts("Sample Call Rate Options");

	subopts.add_options()
		("threshold", po::value<double>(&_thresh),"Drop all markers with a call rate below this threshold")
		;

	opts.add(subopts);
	return opts;
}

void MarkerMissingFilter::parseOptions(const po::variables_map& vm){
	// no need to do anything here!
}

}
}
