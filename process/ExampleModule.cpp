#include "ExampleModule.h" //////CHANGE TO REAL MODULE NAME

#include "util/DataSet.h"

#include <iostream>

using std::cout;
using std::endl;
using std::string;

using Methods::DataSet;

namespace po=boost::program_options;

namespace ProcessLib{

const string ExampleModule::stepname = ExampleModule::doRegister("example");

//PrintSummary()
//used to output results after processing the data
void ExampleModule::PrintSummary(){
	//output for process goes here
	//
	cout << "Hi, I'm in the ExampleModule PrintSummary method\n";
}

//process()
//main method to get the process going and to the work
void ExampleModule::process(DataSet& ds){
	//main processing of step goes here.
	cout << "Hi, I'm in the ExampleModule process method\n";
	cout << "I'm going to process " << ds.num_loci() << " markers, " << ds.num_pedigrees() << " families, and " << ds.num_inds() << " samples!\n";

}

po::options_description& ExampleModule::appendOptions(po::options_description& opts){
	po::options_description subopts("Example Options");

	subopts.add_options()
		("arg", po::value<string>(&arg_string),"Example string argument")
		("switch", po::bool_switch(&arg_bool), "Example switch");

	opts.add(subopts);
	return opts;
}

void ExampleModule::parseOptions(const po::variables_map& vm){
	cout << "Parsing the options for the example module" << endl;
	cout << "NOTE: all variables should be set already, this is just meant "
		 << "to do sophisticated post-processing beyond setting variables\n";
	cout << "You values are:\n";
	cout << "\targ: " << arg_string << std::endl;
	cout << "\tswitch: " << arg_bool << std::endl;
}

}
