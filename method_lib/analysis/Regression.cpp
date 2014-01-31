/*
 * Regression.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#include "Regression.h"

#include "InputManager.h"

#include "util/Logger.h"
#include "util/Marker.h"

#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstring>
#include <limits>

#include <boost/algorithm/string.hpp>

using std::string;
using std::vector;
using std::set;
using std::map;
using std::numeric_limits;

namespace po=boost::program_options;

using Methods::DataSet;

namespace Methods{

namespace Analysis{

Regression::~Regression() {
	for (unsigned int i = 0; i < results.size(); i++){
		delete results[i];
	}
	results.clear();
	out_f.close();
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
		("models", po::value<vector<string> >(&model_files)->composing(), "List of files containing models to generate")
		("incl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to include")
		("excl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to exclude")
		("output", po::value<string>(&out_fn)->default_value("output.txt"), "Name of the file to output results")
		("encoding", po::value<EncodingModel>(&encoding)->default_value("additive"), "Encoding model to use in the regression (additive, dominant, recessive, categorical)")
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

	out_f.open(out_fn.c_str());
	if(!out_f){
		throw std::logic_error("Cannot open regression output file '" + out_fn + "'");
	}
}

void Regression::runRegression(const DataSet& ds){
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

	set<string>::iterator trait_itr = incl_traits.begin();
	while (trait_itr != incl_traits.end()) {
		if (!ds.isTrait(*trait_itr)) {
			Utility::Logger::log_err("WARNING: '" + *trait_itr
					+ "' is not a recognized trait, ignoring.");
			incl_traits.erase(trait_itr++);
		} else {
			++trait_itr;
		}
	}

	set<string>::iterator covar_itr = covar_names.begin();
	while (covar_itr != covar_names.end()) {
		if (!ds.isTrait(*covar_itr)) {
			Utility::Logger::log_err("WARNING: '" + *covar_itr
					+ "' is not a recognized covariate, ignoring.");
			covar_names.erase(covar_itr++);
		} else {
			++covar_itr;
		}
	}

	if(!ds.isTrait(outcome_name)){
		Utility::Logger::log_err("ERROR: '" + outcome_name + "' is not a recognized trait, aborting.", true);
	}

	string first_model = "";
	if(_models.size() > 0){
		first_model=*(_models.begin());
	}

	// Now, set up the outcome variable and the covariates
	DataSet::const_sample_iterator si = ds.beginSample();
	while(si != ds.endSample()){

		// set up the outcome vector
		if(outcome_name.size() == 0){
			if((*si)->isAffectedKnown()){
				_pheno.push_back((*si)->isAffected());
			}else{
				_pheno.push_back(numeric_limits<float>::quiet_NaN());
			}
		}else{
			_pheno.push_back(ds.getTrait(outcome_name,*si));
		}

		// now, set up the covariate vector(s)
		_covars.resize(covar_names.size());
		set<string>::const_iterator covar_itr = covar_names.begin();
		vector<vector<float> >::iterator val_itr = _covars.begin();
		while(covar_itr != covar_names.end()){
			(*val_itr).push_back(ds.getTrait(*covar_itr, *si));
			++covar_itr;
			++val_itr;
		}
	}

	this->initData(first_model, ds);

	Model* nm;
	ModelGenerator mg(ds, incl_traits, pairwise, exclude_markers);
	while( (nm = mg()) ){
		results.push_back(run(nm, ds));
	}

	printResults();
}

Regression::Model* Regression::parseModelStr(const std::string& model_str, const Methods::DataSet& ds) {
	vector<string> model_elements;
	boost::algorithm::split(model_elements, model_str, boost::is_any_of(" \t"), boost::token_compress_on);

	Model* model = new Model();

	for(unsigned int i=0; i < model_elements.size(); i++){
		Marker* const m = ds.getMarker(model_elements[i]);
		if(m != 0){
			model->markers.push_back(m);
		}else if(ds.isTrait(model_elements[i])){
			model->traits.push_back(model_elements[i]);
		}else{ // unrecognized model element!!
			model->markers.clear();
			model->traits.clear();
			break;
		}
	}

	return model;
}

float Regression::getCategoricalWeight(const Marker* m, const DataSet& ds){

	map<const Marker*, float>::const_iterator w_itr = categ_weight.find(m);
	if(w_itr != categ_weight.end()){
		return (*w_itr).second;
	}

	// OK, at this point, we know we don't have the categorical value

	Model mod;
	mod.markers.push_back(m);
	Result* r = run(&mod, ds, true);

	float toret = (categ_weight[m] = r->coeffs[0] / (r->coeffs[0] + r->coeffs[1]));

	delete r;

	return toret;

}

Regression::Result* Regression::run(Model* m, const DataSet& ds, bool categorical) {

	unsigned int numLoci = m->markers.size();
	unsigned int numCovars = covar_names.size();
	unsigned int numTraits = m->traits.size();
	unsigned int n_vars = numLoci + numTraits;

	// determine size of row for each sample in dataset
	unsigned int n_cols = 1 + numLoci + numTraits + numCovars ;
	if(categorical){
		n_cols += numLoci - numTraits;
	} else if (interactions){
		n_cols += (n_vars * (n_vars - 1))/2;
	}

	// Allocate a huge amount of memory for the regression here
	float regress_data[_pheno.size()][n_cols];
	float row_data[n_cols];

	DataSet::const_sample_iterator si = ds.beginSample();
	unsigned int n_samples = 0;
	unsigned int n_missing = 0;

	while(si != ds.endSample()){
		unsigned int pos=0;
		// 1st col is always the phenotype
		row_data[pos++] = _pheno[n_samples];

		// Now, work through the markers
		for (unsigned int i = 0; i<numLoci; i++){
			unsigned char geno = (*si)->getAdditiveGeno(*(m->markers[i]));
			if(geno == Sample::missing_allele){
				row_data[pos++] = numeric_limits<float>::quiet_NaN();
				if(categorical){
					row_data[pos++] = numeric_limits<float>::quiet_NaN();
				}
			} else {
				if(categorical){
					row_data[pos++] = EncodingModel(Methods::Analysis::DOMINANT)(geno);
					row_data[pos++] = EncodingModel(Methods::Analysis::RECESSIVE)(geno);
				} else if (encoding == Methods::Analysis::CATEGORICAL){
					float w = getCategoricalWeight(m->markers[i], ds);
					row_data[pos++] = 1 + (geno==1) * w + (geno==2);
				} else {
					row_data[pos++] = encoding(geno);
				}
			}
		}

		// Now, work through the traits
		for(unsigned int i = 0; i < numTraits * (!categorical); i++){
			row_data[pos++] = ds.getTrait(m->traits[i], *si);
		}

		// Now, the interaction terms
		for(unsigned int i=0; i < (n_vars) * interactions * (!categorical); i++){
			for(unsigned int j=i+1; j < (n_vars); j++){
				row_data[pos++] = row_data[i+1]*row_data[j+1];
			}
		}

		// Now, the covariates
		for(unsigned int i=0; i<_covars.size(); i++){
			row_data[pos++] = _covars[i][n_samples];
		}

		// OK, check for missingness before adding it to the dataset
		bool ismissing = false;
		for(unsigned int i=0; i<pos; i++){
			ismissing |= std::isnan(row_data[i]);
		}

		if(!ismissing){
			// do a fast memory copy into the data for regression
			std::memcpy(regress_data[n_samples - n_missing], row_data, sizeof(float) * (pos));
		} else {
			++n_missing;
		}

		++n_samples;
	}

	delete m;

	return calculate(regress_data[0], n_cols, n_samples - n_missing);

}

Regression::Model* Regression::ModelGenerator::operator()() {

	Model* m = 0;
	if (_targeted) {
		// models are given one per line here
		if(_mitr != _mend){
			m = Regression::parseModelStr(*_mitr, _ds);
			++_mitr;
		}
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

