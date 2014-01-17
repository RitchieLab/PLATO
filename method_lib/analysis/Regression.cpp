/*
 * Regression.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#include "Regression.h"

#include "InputManager.h"

using std::string;
using std::vector;

namespace po=boost::program_options;


namespace Methods{

namespace Analysis{

po::options_description& Regression::addOptions(po::options_description& opts){
	po::options_description regress_opts("Regression Options");

	regress_opts.add_options()
		("covariates", po::value<vector<string> >()->composing(),
			"A list of covariates to use in the model")
		("outcome", po::value<string>(&outcome_name)->default_value(""),
			"Use a given covariate as the regression model outcome")
		("exclude-markers", po::bool_switch(&exclude_markers),
			"Do not include markers in the generated models (use for EWAS)")
		("use-traits", po::bool_switch(&include_traits),
			"Include the traits not used as covariates in the generated models")
		("interactions", po::bool_switch(&interactions),
			"Include interactions in the models generated")
		("pairwise", po::bool_switch(&pairwise),
			"When auto-generating models, include two variables exhaustively")
		("correction", po::value<vector<string> >()->composing(), "p-value correction method(s)");

	opts.add(regress_opts);

	return opts;
}

void Regression::parseOptions(const boost::program_options::variables_map& vm){
	if(vm.count("correction")){
		InputManager::parseInput(vm["correction"].as<vector<string> >(), corr_methods);
	}

	if(vm.count("covariates")){
		InputManager::parseInput(vm["covariates"].as<vector<string> >(), covar_names);
	}
}


Regression::Regression() {
	// TODO Auto-generated constructor stub

}

Regression::~Regression() {
	// TODO Auto-generated destructor stub
}

}

}

