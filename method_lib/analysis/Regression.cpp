/*
 * Regression.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#include "Regression.h"

#include "util/InputManager.h"
#include "util/Logger.h"

#include "data/Marker.h"
#include "data/Sample.h"

#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <iostream>

#include <gsl/gsl_cdf.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using std::string;
using std::vector;
using std::set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::min;
using std::max;

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Data::Sample;
using PLATO::Utility::InputManager;

namespace po=boost::program_options;


namespace PLATO{

namespace Analysis{

Regression::~Regression() {
	for (unsigned int i = 0; i < results.size(); i++){
		delete results[i];
	}
	results.clear();

	map<const Marker*, Result*>::const_iterator mu_itr = _marker_uni_result.begin();
	while(mu_itr != _marker_uni_result.end()){
		delete (*mu_itr).second;
		++mu_itr;
	}
	_marker_uni_result.clear();

	map<string, Result*>::const_iterator tu_itr = _trait_uni_result.begin();
	while(tu_itr != _trait_uni_result.end()){
		delete (*tu_itr).second;
		++tu_itr;
	}
	_trait_uni_result.clear();

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
		("seaparator", po::value<string>(&sep)->default_value("\t", "<TAB>"), "Separator to use when outputting results file")
		("encoding", po::value<EncodingModel>(&encoding)->default_value("additive"), "Encoding model to use in the regression (additive, dominant, recessive, categorical)")
		("show-univariate", po::bool_switch(&show_uni), "Show univariate results in multivariate models")
		;

	opts.add(regress_opts);

	return opts;
}

void Regression::printVarHeader(const string& var_name){
	out_f << var_name << "_Pval" << sep
		  << var_name << "_beta" << sep
		  << var_name << "_SE" << sep;
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

	string model_str = "";
	if(_models.size() > 0){
		model_str=*(_models.begin());
	}

	unsigned int n_snp = 0;
	unsigned int n_trait = 0;

	if(model_str.size() > 0){
		Model* m = Regression::parseModelStr(model_str, ds);
		n_snp = m->markers.size();
		n_trait = m->markers.size();
		delete m;
	} else{
		// Calculate the number of SNPs / Env vars based on the options passed
		// in to the Regression object.
		n_snp = (!exclude_markers) * (1 + pairwise);
		n_trait = (incl_traits.size() > 0) * (1 + exclude_markers);
	}

	// Now, print the header
	for(unsigned int i=0; i<n_snp; i++){
		out_f << "Var" << i+1 << "_ID" << sep
			  << "Var" << i+1 << "_Pos" << sep
			  << "Var" << i+1 << "_MAF" << sep;
	}

	for(unsigned int i=0; i<n_trait; i++){
		out_f << "Var" << n_snp + i + 1 << "_ID" << sep;
	}

	out_f << "N_Missing" << sep;

	if(n_snp + n_trait > 1){
		for(unsigned int i=0; show_uni && i < n_snp + n_trait; i++){
			printVarHeader("Uni_Var" + boost::lexical_cast<string>(i+1));
		}

		if(interactions){
			for(unsigned int i=0; i < n_snp + n_trait; i++){
				printVarHeader("Red_Var"+ boost::lexical_cast<string>(i+1));
			}

			for (unsigned int i = 0; i < n_snp + n_trait; i++) {
				printVarHeader("Full_Var" + boost::lexical_cast<string>(i+1));
			}

			for(unsigned int i=0; i<n_snp+n_trait; i++){
				for(unsigned int j=1; j<n_snp+n_trait; j++){
					printVarHeader("Full_Var" + boost::lexical_cast<string>(i+1)
							+ "_Var" + boost::lexical_cast<string>(j+1));
				}
			}

			out_f << "Red_Model_Pval" << sep << "Full_Model_Pval" << sep;
		}
	}

	for(unsigned int i=0; !interactions && i < n_snp + n_trait; i++){
		printVarHeader("Var" + boost::lexical_cast<string>(i+1));
	}

	// If we are looking at single SNP models w/ categorical weight, print said weight!
	if(encoding == Encoding::CATEGORICAL && n_snp == 1 && n_trait == 0){
		out_f << "Categ_Weight" << sep;
	}

	out_f << "Overall_Pval";

	set<CorrectionModel>::const_iterator c_itr = corr_methods.begin();
	while(c_itr != corr_methods.end()){
		out_f << sep << "Overall_Pval_adj_" << *c_itr;
	}

	out_f << std::endl;

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

	this->initData(model_str, ds);

	Model* nm;
	ModelGenerator mg(ds, incl_traits, pairwise, exclude_markers);
	while( (nm = mg()) ){
		Result* r = run(nm, ds, false);

		// Run some univariate models
		if(show_uni && n_snp + n_trait > 1){
			// Here, we are going to reverse the order of marker/trait variables
			// because we're adding the results onto the FRONT of the model

			for(unsigned int i=nm->traits.size(); i>0; i--){
				map<string, Result*>::const_iterator t_itr = _trait_uni_result.find(nm->traits[i-1]);
				// If this is the case, we have not seen this marker before
				if(t_itr == _trait_uni_result.end()){
					Model m;
					m.traits.push_back(nm->traits[i-1]);
					t_itr = _trait_uni_result.insert(_trait_uni_result.begin(),
							std::make_pair(nm->traits[i-1], run(&m, ds)));
				}

				// Add the result to the front of the deque
				r->coeffs.push_front((*t_itr).second->coeffs[0]);
				r->p_vals.push_front((*t_itr).second->p_vals[0]);
				r->stderr.push_front((*t_itr).second->stderr[0]);
			}

			for(unsigned int i=nm->markers.size(); i>0; i--){
				map<const Marker*, Result*>::const_iterator m_itr = _marker_uni_result.find(nm->markers[i-1]);
				// If this is the case, we have not seen this marker before
				if(m_itr == _marker_uni_result.end()){
					Model m;
					m.markers.push_back(nm->markers[i-1]);
					m_itr = _marker_uni_result.insert(_marker_uni_result.begin(),
							std::make_pair(nm->markers[i-1], run(&m, ds)));
				}

				// Add the result to the front of the deque
				r->coeffs.push_front((*m_itr).second->coeffs[0]);
				r->p_vals.push_front((*m_itr).second->p_vals[0]);
				r->stderr.push_front((*m_itr).second->stderr[0]);
			}
		}

		if(interactions){
			Result* r_full = run(nm, ds, true);

			stringstream ss;
			ss << r->p_val << sep << r_full->p_val << sep;
			for(unsigned int i=0; i<r_full->coeffs.size(); i++){
				r->coeffs.push_back(r_full->coeffs[i]);
				r->stderr.push_back(r_full->stderr[i]);
				r->p_vals.push_back(r_full->p_vals[i]);
			}

			// calculate a new log likelihood
			r->log_likelihood = -2 * (r->log_likelihood - r_full->log_likelihood);
			r->p_val = gsl_cdf_chisq_Q(r->log_likelihood,1);
			delete r_full;

		}

		results.push_back(r);
		delete nm;
	}

	printResults();
}

Regression::Model* Regression::parseModelStr(const std::string& model_str, const DataSet& ds) {
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
	Result* r = run(&mod, ds, false, true);

	// NOTE: this assigns to the map at the same step
	float toret = (categ_weight[m] = r->coeffs[0] / (r->coeffs[0] + r->coeffs[1]));

	delete r;

	return toret;

}

Regression::Result* Regression::run(const Model* m, const DataSet& ds, bool interact, bool categorical) {

	unsigned int numLoci = m->markers.size();
	unsigned int numCovars = covar_names.size();
	unsigned int numTraits = m->traits.size();
	unsigned int n_vars = numLoci + numTraits;
	unsigned int n_kept = n_vars;

	// determine size of row for each sample in dataset
	unsigned int n_cols = 1 + numLoci + numTraits + numCovars ;
	if(categorical){
		n_cols += numLoci - numTraits;
	} else if (interactions){
		n_kept += (n_vars * (n_vars - 1))/2;
		n_cols += n_kept - n_vars;
	}

	// Allocate a huge amount of memory for the regression here
	double regress_data[_pheno.size()][n_cols];
	double row_data[n_cols];
	unsigned char geno[numLoci];
	float maf_sum[numLoci];

	DataSet::const_sample_iterator si = ds.beginSample();
	unsigned int n_samples = 0;
	unsigned int n_missing = 0;

	double categ_weight[numLoci];
	if (encoding == Encoding::CATEGORICAL){
		for (unsigned int i = 0; i<numLoci; i++){
			categ_weight[i] = getCategoricalWeight(m->markers[i], ds);
		}
	}

	while(si != ds.endSample()){
		unsigned int pos=0;
		// 1st col is always the phenotype
		row_data[pos++] = _pheno[n_samples];

		// Now, work through the markers
		for (unsigned int i = 0; i<numLoci; i++){
			geno[i] = (*si)->getAdditiveGeno(*(m->markers[i]));
			if(geno[i] == Sample::missing_allele){
				row_data[pos++] = numeric_limits<double>::quiet_NaN();
				if(categorical){
					row_data[pos++] = numeric_limits<double>::quiet_NaN();
				}
			} else {
				if(categorical){
					row_data[pos++] = EncodingModel(Encoding::DOMINANT)(geno[i]);
					row_data[pos++] = EncodingModel(Encoding::RECESSIVE)(geno[i]);
				} else if (encoding == Encoding::CATEGORICAL){
					// Returns:
					// {0,w,1}    , w in [0,1]
					// {1,1-w,0}  , w in [-1,0]
					// {0,1,1/w}  , w in [1, inf)
					// {1,0,1-1/w}, w in [-inf, -1]
					row_data[pos++] = (categ_weight[i] < 0)
							        + max(min(1.0, categ_weight[i]),-1.0) * (geno[i]==1)
							        + max(min(1.0, 1/categ_weight[i]),-1.0) * (geno[i]==2);
				} else {
					row_data[pos++] = encoding(geno[i]);
				}
			}
		}

		// Now, work through the traits
		for(unsigned int i = 0; i < numTraits * (!categorical); i++){
			row_data[pos++] = ds.getTrait(m->traits[i], *si);
		}

		// Now, the interaction terms
		for(unsigned int i=0; i < (n_vars) * interact * (!categorical); i++){
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
			for(unsigned int i=0; i<numLoci; i++){
				maf_sum[i] += geno[i];
			}
			// do a fast memory copy into the data for regression
			std::memcpy(regress_data[n_samples - n_missing], row_data, sizeof(double) * (pos));
		} else {
			++n_missing;
		}

		++n_samples;
	}

	Result* r = calculate(regress_data[0], n_cols, n_samples - n_missing);
	r->coeffs.resize(n_kept);
	r->stderr.resize(n_kept);
	r->p_vals.resize(n_kept);

	stringstream ss;

	for(unsigned int i=0; i<numLoci; i++){
		ss << m->markers[i]->getID() << sep
		   << m->markers[i]->getChromStr() << ":" << m->markers[i]->getLoc()
		   << sep << maf_sum[i] / (2*static_cast<float>(n_samples-n_missing)) << sep;
	}
	for(unsigned int i=0; i<numTraits; i++){
		ss << m->traits[i] << sep;
	}

	// print the # missing from this model
	ss << n_missing << sep;

	r->prefix = ss.str();

	ss.clear();
	if(encoding == Encoding::CATEGORICAL && m->markers.size() == 1 && m->traits.size() == 0 && !categorical){
		ss << categ_weight << sep;
		r->suffix += ss.str();
	}

	return r;
}

void Regression::printResults(){
	// Perhaps sort the deque of results based on overall p-values

	std::sort(results.begin(), results.end(), result_sorter());

	// Now, let's do some multiple test correction!
	map<CorrectionModel, vector<float> > pval_corr;

	// don't bother with the work if there's no correction desired!
	if (corr_methods.size()) {

		vector<float> pv_in;
		pv_in.reserve(results.size());
		for (unsigned int i = 0; i < results.size(); i++) {
			pv_in.push_back(results[i]->p_val);
		}

		set<CorrectionModel>::const_iterator c_itr = corr_methods.begin();
		while(c_itr != corr_methods.end()){
			Correction::getCorrectionMethod(*c_itr)->correct(pv_in, pval_corr[(*c_itr)]);
		}

	}


	for(unsigned int i=0; i<results.size(); i++){
		// If we've gone over the threshold, stop printing!
		if(results[i]->p_val > cutoff_p){
			break;
		}
		// print a single line
		out_f << results[i]->prefix;

		for(unsigned int j=0;j<results[i]->coeffs.size(); j++){
			out_f << results[i]->p_vals[j] << sep
				  << results[i]->coeffs[j] << sep
				  << results[i]->stderr[j] << sep;
		}

		out_f << results[i]->suffix << results[i]->p_val;

		// print some creected p-values!
		// NOTE: we assume that "set" and "map" use the same ordering!!
		map<CorrectionModel, vector<float> >::const_iterator p_itr = pval_corr.begin();
		while(p_itr != pval_corr.end()){
			out_f << sep << (*p_itr).second[i];
		}

		out_f << std::endl;
	}
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

