#include "RecodeAlleles.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"

#include "util/Logger.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

#include <boost/algorithm/string.hpp>

using std::ifstream;
using std::ofstream;
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

po::options_description& RecodeAlleles::appendOptions(po::options_description& opts){
	po::options_description subopts("Example Options");

	subopts.add_options()
		("file", po::value<string>(&recode_fn),"File of markers and given referent allele")
		("input-map", po::bool_switch(&map_input), "Input file of markers is in (extended) map format")
		("map3", po::bool_switch(&map3), "Use 3-column map format for input and output")
		("auto", po::bool_switch(&auto_recode), "Automatically make the referent allele the major allele")
		("out", po::value<string>(&out_fn), "Output a file of markers along with their referent allele")
		("output-map", po::bool_switch(&map_output), "Print output in (extended) map format")
		;

	opts.add(subopts);
	return opts;
}

void RecodeAlleles::parseOptions(const po::variables_map& vm){
	if(auto_recode && recode_fn.size() > 0){
		Logger::log_err("WARNING: --file and --auto options both specified: --file will take precedence");
	}

	if(map_input && recode_fn.size() == 0){
		Logger::log_err("WARNING: --input-map given without a --file option: ignoring --input-map option");
		map_input = false;
	}

	if(map_output && out_fn.size() == 0){
		Logger::log_err("WARNING: --output-map given without a --out option: ignoring --output-map");
	}

	if(map3 && !(map_input || map_output)){
		Logger::log_err("WARNING: --map3 given without --map-input or --map-output: ignoring --map3");
	}
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
				string chr, pos, id, ref;
				if(map_input){
					ss >> chr;
				}
				ss >> id;
				if(map_input){
					ss >> pos;
					if(!map3){
						// I'm just going to overwrite this, so no worries
						// about creating a temp. variable here
						ss >> ref;
					}
				}
				ss >> ref;
				Marker* m = ds.getMarker(id);
				if(!m && map_input){
					m = ds.getMarker(chr, atoi(pos.c_str()));
				}
				if(!m){
					Logger::log_err("WARNING: could not find marker with ID: '" + id + "'");
				} else if (!m->setRefAllele(ref)) {
					Logger::log_err("WARNING: could not set reference allele '" + ref + "' for marker '" + id + "'");
				}
			}
		}

		rf.close();
	}

	if(out_fn != ""){
		ofstream outf(out_fn.c_str());
		DataSet::const_marker_iterator cmi = ds.beginMarker();
		while(cmi != ds.endMarker()){

			if(map_output){
				outf << (*cmi)->getChromStr() << "\t";
			}
			outf << (*cmi)->getID() << "\t";
			if(map_output){
				outf << (*cmi)->getLoc() << "\t";
				if(!map3){
					outf << "0\t";
				}
			}
			outf << (*cmi)->getRefAllele() << std::endl;

		}

		outf.close();
	}
}

}
}
