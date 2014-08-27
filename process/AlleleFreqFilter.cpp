#include "AlleleFreqFilter.h" //////CHANGE TO REAL MODULE NAME

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

const string AlleleFreqFilter::stepname = AlleleFreqFilter::doRegister("filter-maf");

//PrintSummary()
//used to output results after processing the data
void AlleleFreqFilter::PrintSummary(){
	//output for process goes here
	//
	Utility::Logger::log("Allele Frequency Filter: Dropped " + boost::lexical_cast<string>(_n_filtered) + " Markers.\n");
}

//process()
//main method to get the process going and to the work
void AlleleFreqFilter::process(DataSet& ds){

	_n_filtered = 0;
	DataSet::marker_iterator mi = ds.beginMarker();
	float curr_maf = 0;
	while(mi != ds.endMarker()){
		curr_maf = (*mi)->issetMAF() ? (*mi)->getMAF() : (*mi)->calcMAF(ds);
		if(curr_maf < _min_thresh || curr_maf > _max_thresh){
			(*mi)->setEnabled(false);
			++_n_filtered;
		}
		++mi;
	}
}

po::options_description& AlleleFreqFilter::appendOptions(po::options_description& opts){
	po::options_description subopts("Sample Call Rate Options");

	subopts.add_options()
		("min", po::value<double>(&_min_thresh)->default_value(0.05, "0.05"),"Drop markers with MAF below this")
		("max", po::value<double>(&_max_thresh)->default_value(1, "1.0"),"Drop markers with MAF above this")
		;

	opts.add(subopts);
	return opts;
}

void AlleleFreqFilter::parseOptions(const po::variables_map& vm){
	// no need to do anything here!
}

}
}
