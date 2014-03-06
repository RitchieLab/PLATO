#include "ExampleModule.h" //////CHANGE TO REAL MODULE NAME

#include "data/DataSet.h"

#include <iostream>

#include <gsl/gsl_multifit.h>

using std::cout;
using std::endl;
using std::string;

using PLATO::Data::DataSet;

namespace po=boost::program_options;

namespace PLATO{
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

	// Test some GSL regression code here - easier to debug in Eclipse
	gsl_matrix* dat = gsl_matrix_calloc(6,3);
	gsl_vector* val = gsl_vector_calloc(6);
	gsl_matrix_set_identity(dat);
	// make the 2nd column all 0s (singularity)
	gsl_matrix_set(dat,1,1,0);
	// Make the first column all 1
	for(int i=0; i<6; i++){
		gsl_matrix_set(dat,i,0,1);
	}

	gsl_vector_set_basis(val,0);
	gsl_vector_set(val,1,1);
	gsl_vector_set(val,2,1);

	// Now we have:
	// val = [1 1 1 0 0 0]
	// mat = [1 1 1 1 1 1; 0 0 0 0 0 0; 0 0 1 0 0 0]'

	// try to run a linear regression:
	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc (6, 3);

	gsl_vector* c = gsl_vector_alloc(3);
	gsl_matrix* cov = gsl_matrix_alloc(3,3);
	double chisq;

	gsl_multifit_linear(dat, val, c, cov, &chisq, ws);

	gsl_matrix_free(cov);
	gsl_vector_free(c);
	gsl_multifit_linear_free(ws);
	gsl_matrix_free(dat);
	gsl_vector_free(val);

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
}
