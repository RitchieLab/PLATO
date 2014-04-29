#include "OutputBeagle.h" //////CHANGE TO REAL MODULE NAME

#include "util/Logger.h"

#include <fstream>

using std::string;
using std::ofstream;
using std::endl;

using PLATO::Data::DataSet;
using PLATO::Data::Family;
using PLATO::Data::Sample;
using PLATO::Utility::Logger;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string OutputBeagle::stepname = OutputBeagle::doRegister("output-beagle");

po::options_description& OutputBeagle::appendOptions(po::options_description& opts){
	po::options_description subopts("BEAGLE Output Options");

	subopts.add_options()
		("prefix", po::value<string>(&_prefix)->default_value("output"), "Prefix for output data")
		("suffix", po::value<string>(&_suffix)->default_value("bgl"), "Suffix for genotype data")
		("marker-suffix", po::value<string>(&_marker_suff)->default_value("markers"), "Suffix for marker data")
		("incl-traits", po::bool_switch(&_incl_trait), "Include traits in the genotype output data")
		("pair", po::bool_switch(&_pair), "Data to output is pair data")
		("trio", po::bool_switch(&_trio), "Data to output is trio data")
		;

	opts.add(subopts);
	return opts;
}

void OutputBeagle::parseOptions(const po::variables_map& vm){

	if(_pair && _trio){
		Logger::log_err("ERROR: Data cannot be both pair and trio data", true);
	}

}


//process()
//main method to get the process going and to the work
void OutputBeagle::process(DataSet& ds){
	//main processing of step goes here.

	// make sure the markers are sorted
	ds.sortMarkers();
	string prevChr = "";

	DataSet::const_marker_iterator mi = ds.beginMarker();
	SampleGenerator* sg = getSampleGenerator(ds);
	const Sample* s = 0;
	ofstream genof;
	ofstream markerf;

	while(mi != ds.endMarker()){
		if((*mi)->getChromStr() != prevChr){
			prevChr = (*mi)->getChromStr();
			if(genof){
				genof.close();
			}
			if(markerf){
				markerf.close();
			}
			genof.open((_prefix + "." + prevChr + "." + _suffix).c_str());
			markerf.open((_prefix + "." + prevChr + "." + _marker_suff).c_str());

			genof << "I" << _sep << "id";

			sg->reset();
			while( (s = sg->next()) != 0){
				genof << _sep << s->getID() << _sep << s->getID();
			}
			genof << endl;

			if(_incl_trait){
				genof << "A" << _sep << "status";
				sg->reset();
				while( (s = sg->next()) != 0){
					for(unsigned int i=0; i<2; i++){
						genof << _sep << s->isAffectedKnown() * (s->isAffected() + 1);
					}
				}
				genof << endl;

				DataSet::const_trait_iterator ti = ds.beginTrait();
				while(ti != ds.endTrait()){
					genof << "T" << _sep << *ti;
					sg->reset();
					while( (s = sg->next()) != 0){
						for(unsigned int i=0; i<2; i++){
							genof << _sep << ds.getTrait(*ti, s);
						}
					}
					genof << endl;
					++ti;
				}
			}
		}

		markerf << (*mi)->getID() << _sep << (*mi)->getLoc();
		for(unsigned int i=0; i<(*mi)->getNumAlleles(); i++){
			markerf << _sep << (*mi)->getAllele(i);
		}
		markerf << endl;

		genof << "M" << _sep << (*mi)->getID();

		sg->reset();
		while( (s = sg->next()) != 0){
			std::pair<unsigned char, unsigned char> geno = s->getGeno(**mi);
			genof << _sep << ((*mi)->isAlleleMissing(geno.first) ?
						string(".") : (*mi)->getAllele(geno.first))
				  << _sep << ((*mi)->isAlleleMissing(geno.second) ?
						string(".") : (*mi)->getAllele(geno.second));
		}
		genof << endl;

		++mi;
	}

	genof.close();
	markerf.close();
	delete sg;

}

OutputBeagle::SampleGenerator* OutputBeagle::getSampleGenerator(const DataSet& ds){
	if(_trio || _pair){
		return new FamilySampleGenerator(ds, 1 + _trio);
	} else {
		return new BasicSampleGenerator(ds);
	}
}

OutputBeagle::FamilySampleGenerator::FamilySampleGenerator(const DataSet& ds, unsigned int n_found)
	: OutputBeagle::SampleGenerator(ds), _FOUNDERS(n_found), _curr_founder(0), _fi(ds.beginFamily()) {

	// let's check to make sure everything is as it should be!
	while(_fi != ds.endFamily()){
		Family::const_member_iterator fmi = (*_fi)->beginMembers(false);

		_curr_founder = 0;
		while(_curr_founder < _FOUNDERS){
			++_curr_founder;
			if(fmi == (*_fi)->endMembers()){
				Logger::log_err("ERROR: Not enough founders found in family '"+ (*_fi)->getFamID() + "'", true);
			}
			++fmi;
		}

		fmi = (*_fi)->beginMembers(false);
		Sample::const_child_iterator ci = (*fmi)->beginChild();
		unsigned int n_kids = 0;
		while(ci != (*fmi)->endChild()){
			++n_kids;
			++ci;
		}
		if(n_kids != 1){
			Logger::log_err("WARNING: Multiple children found in family '" + (*_fi)->getFamID() + "', taking only the first child.");
		}

		++_fi;
	}

	_fi = ds.beginFamily();
	_curr_founder = 0;

}

const Sample* OutputBeagle::FamilySampleGenerator::next(){
	const Sample* toret = 0;

	if(_curr_founder == _FOUNDERS){
		toret = *((*((*_fi)->beginMembers(false)))->beginChild());
		_curr_founder = 0;
		++_fi;
	} else if(_fi == _ds.endFamily()) {
		Family::const_member_iterator fmi = (*_fi)->beginMembers(false);
		for(unsigned int i=0; i<_curr_founder; i++){
			++fmi;
		}
		toret = (*fmi);
		++_curr_founder;
	}

	return toret;
}


}
}
