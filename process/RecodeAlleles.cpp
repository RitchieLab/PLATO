#include "RecodeAlleles.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"

#include "util/Logger.h"

#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>

using std::ifstream;
using std::stringstream;
using std::string;

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Utility::Logger;


namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string RecodeAlleles::stepname = RecodeAlleles::doRegister("recode-alleles");

//PrintSummary()
//used to output results after processing the data
void RecodeAlleles::PrintSummary(){
	//output for process goes here
}

//process()
//main method to get the process going and to the work
void RecodeAlleles::process(DataSet& ds){
	if(auto_recode){
		DataSet::marker_iterator mi = ds.beginMarker();
		while(mi != ds.endMarker()){
			(*mi)->setRefAlleleIdx((*mi)->majorAlleleIdx(ds));
			++mi;
		}
	}
	if(recode_fn != ""){
		ifstream rf(recode_fn.c_str());
		string l;
		while(getline(rf, l)){
			boost::algorithm::trim(l);
			if(l.size() > 0 && l[0] != '#'){
				stringstream ss(l);
				string id, ref;
				ss >> id >> ref;
				Marker* m = ds.getMarker(id);
				if(!m){
					Logger::log_err("WARNING: could not find marker with ID: '" + id + "'");
				} else if (!m->setRefAllele(ref)) {
					Logger::log_err("WARNING: could not set reference allele '" + ref + "' for marker '" + id + "'");
				}
			}
		}
	}
}

po::options_description& RecodeAlleles::appendOptions(po::options_description& opts){
	po::options_description subopts("Example Options");

	subopts.add_options()
		("file", po::value<string>(&recode_fn),"File of markers and given referent allele")
		("auto", po::bool_switch(&auto_recode), "Automatically make the referent allele the major allele");

	opts.add(subopts);
	return opts;
}

void RecodeAlleles::parseOptions(const po::variables_map& vm){
	if(auto_recode && vm.count("file")){
		Logger::log_err("WARNING: --file and --auto options both specified: --file will take precedence");
	}
}

}
}
