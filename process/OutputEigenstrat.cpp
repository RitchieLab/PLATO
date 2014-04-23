#include "OutputEigenstrat.h"

#include "data/DataSet.h"

#include <fstream>

using std::string;
using std::ofstream;

using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const string OutputEigenstrat::stepname = OutputEigenstrat::doRegister("output-eigenstrat");

//PrintSummary()
//used to output results after processing the data
void OutputEigenstrat::PrintSummary(){
	//output for process goes here
	//

}

//process()
//main method to get the process going and to the work
void OutputEigenstrat::process(DataSet& ds){
	//main processing of step goes here.

	const static string sep = "\t";

	// print out the marker (snp) file
	ofstream snp_f(_snp_fn.c_str());
	DataSet::const_marker_iterator mi = ds.beginMarker();
	DataSet::const_marker_iterator me = ds.endMarker();
	while(mi != me){
		snp_f << (*mi)->getID() << sep << (*mi)->getChrom() << sep
			  << 0 << sep << (*mi)->getLoc() << std::endl;
		++mi;
	}
	snp_f.close();

	// print out the individual file
	ofstream indiv_f(_indiv_fn.c_str());
	DataSet::const_sample_iterator si = ds.beginSample();
	DataSet::const_sample_iterator se = ds.endSample();
	while(si != se){
		indiv_f << (*si)->getID() << sep
				<< ((*si)->isGenderKnown() ? ((*si)->isMale() ? "M" : "F") : "U") << sep
				<< ((*si)->isAffectedKnown() * ((*si)->isAffected() + 1)) << std::endl;
		++si;
	}
	indiv_f.close();

	// print out the genotype file
	ofstream geno_f(_geno_fn.c_str());
	mi = ds.beginMarker();
	while(mi != me){
		si = ds.beginSample();
		while(si != se){
			geno_f << ((*si)->isMissing(**mi) ? missing_val : (*si)->getAdditiveGeno(**mi));
			++si;
		}
		geno_f << std::endl;
		++mi;
	}
	geno_f.close();

}

po::options_description& OutputEigenstrat::appendOptions(po::options_description& opts){
	po::options_description subopts("Eigenstrat Output Options");

	subopts.add_options()
		("prefix", po::value<string>(&_prefix)->default_value("plato"), "Prefix for output filenames")
		("genotype", po::value<string>(&_geno_fn),"Genotype file name")
		("snp", po::value<string>(&_snp_fn),"SNP file name")
		("indiv", po::value<string>(&_indiv_fn),"Individual file name")
		;

	opts.add(subopts);
	return opts;
}

void OutputEigenstrat::parseOptions(const po::variables_map& vm){
	if(!vm.count("genotype")){
		_geno_fn = _prefix + ".geno";
	}

	if(!vm.count("snp")){
		_snp_fn = _prefix + ".snp";
	}

	if(!vm.count("indiv")){
		_indiv_fn = _prefix + ".indiv";
	}
}

}
}
