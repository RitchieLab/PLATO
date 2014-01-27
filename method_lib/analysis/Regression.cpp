/*
 * Regression.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#include "Regression.h"

#include "InputManager.h"

#include "util/Marker.h"

#include <algorithm>

using std::string;
using std::vector;
using std::set;

namespace po=boost::program_options;

using Methods::DataSet;

namespace Methods{

namespace Analysis{

Regression::~Regression() {
	for (unsigned int i = 0; i < results.size(); i++){
		delete results[i];
	}
	results.clear();
}

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
		("correction", po::value<vector<string> >()->composing(), "p-value correction method(s)")
		("models", po::value<vector<string> >(&model_files)->composing(), "List of models to generate")
		("incl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to include")
		("excl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to exclude")
		;

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

	if(vm.count("incl-traits")){
		InputManager::parseInput(vm["incl-traits"].as<vector<string> >(), incl_traits);
	}

	if(vm.count("incl-traits")){
		InputManager::parseInput(vm["excl-traits"].as<vector<string> >(), excl_traits);
	}

}

void Regression::runRegression(DataSet& ds){
	if(include_traits){
	// First let's get all of the traits to include...
		if(incl_traits.size() == 0 && include_traits){
			incl_traits.insert(ds.beginTrait(), ds.endTrait());
		}
		set<string> tmp_trait;
		set<string>::iterator tmp_itr = tmp_trait.begin();
		// Now, take out the excluded traits
		std::set_difference(incl_traits.begin(), incl_traits.end(),
							excl_traits.begin(), excl_traits.end(),
							std::inserter(tmp_trait, tmp_trait.begin()));

		incl_traits.clear();

		// Now, take out the covariates
		std::set_difference(tmp_trait.begin(), tmp_trait.end(),
						    covar_names.begin(), covar_names.end(),
						    std::inserter(incl_traits, incl_traits.begin()));
	}else{
		incl_traits.clear();
	}


	Model* nm;
	ModelGenerator mg(ds, incl_traits, pairwise, exclude_markers);
	while( (nm = mg()) ){
		results.push_back(run(nm));
	}
}

Regression::Model* Regression::ModelGenerator::operator()() {

	Model* m = 0;
	if (_targeted) {
		// models are given one per line here


	} else {
		// We must want exhaustive pairwise models
		if (_pairwise) {

			// We want Env vars
			if (_traits) {
				// we want exhaustive EnvxEnv
				if (_nomarker) {
					if (++_ti2 == _tend && ++_titr != _tend) {
						_ti2 = _titr;
						++_ti2;
					}
					if (_titr != _tend && _ti2 != _tend) {
						m = new Model();
						m->traits.push_back(*_titr);
						m->traits.push_back(*_ti2);
					}

				// We want exhaustive SNPxSNPxEnv models
				} else {

					if (_mi2 == _ds.endMarker()) {
						if (++_mi1 != _ds.endMarker()) {
							_mi2 = _mi1;
							if (++_mi2 == _ds.endMarker()) {
								_mi1 = _ds.beginMarker();
								_mi2 = _mi1;
								++_titr;
							}
						}
					}
					if (_titr != _tend) {
						m = new Model();
						m->markers.push_back(*_mi1);
						m->markers.push_back(*_mi2);
						m->traits.push_back(*_titr);
						++_mi2;
					}
				}
			// OK, we just want SNPxSNP models
			} else{

				if (_mi2 == _ds.endMarker()) {
					if (++_mi1 != _ds.endMarker()) {
						_mi2 = _mi1;
						++_mi2;
					}
				}
				if(_mi2 != _ds.endMarker()){
					m = new Model();
					m->markers.push_back(*_mi1);
					m->markers.push_back(*_mi2);
					++_mi2;
				}

			}
		// We want Marker x Trait models (i.e. GxE)
		}else if(_traits){ // Note: !_pairwise == true here
			if(_mi1 == _ds.endMarker()){
				_mi1 = _ds.beginMarker();
				++_titr;
			}
			if(_mi1 != _ds.endMarker() && _titr != _tend) {
				m = new Model();
				m->markers.push_back(*_mi1);
				m->traits.push_back(*_titr);
				++_mi1;
			}
		// Must want single variable models
		}else{
			// Env only
			if(_nomarker){
				if(_titr != _tend){
					m = new Model();
					m->traits.push_back(*_titr);
					++_titr;
				}
			// SNP only
			} else {
				if(_mi1 != _ds.endMarker()){
					m = new Model();
					m->markers.push_back(*_mi1);
					++_mi1;
				}
			}
		}
	}
	// NOTE: m == 0 if no more models can be generated!
	// IN targeted mode, m != 0, but m->markers.size() == 0 && m->traits.size() == 0 in the case of a "bad" model
	return m;

}


}
}

