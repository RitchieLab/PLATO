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
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

using std::string;
using std::vector;
using std::set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::min;
using std::max;
using std::ifstream;

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Data::Sample;
using PLATO::Utility::InputManager;
using PLATO::Utility::Logger;

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
		("correction", po::value<vector<string> >()->composing(), ("p-value correction method(s) (" + Correction::listCorrectionMethods() + ")").c_str())
		("models", po::value<vector<string> >(&model_files)->composing(), "List of files containing models to generate")
		("incl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to include")
		("excl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to exclude")
		("output", po::value<string>(&out_fn)->default_value("output.txt"), "Name of the file to output results")
		("seaparator", po::value<string>(&sep)->default_value("\t", "<TAB>"), "Separator to use when outputting results file")
		("encoding", po::value<EncodingModel>(&encoding)->default_value("additive"), "Encoding model to use in the regression (additive, dominant, recessive, weighted, codominant)")
		("show-univariate", po::bool_switch(&show_uni), "Show univariate results in multivariate models")
		("thresh", po::value<float>(&cutoff_p)->default_value(1.0f), "Threshold for printing resultant models")
		("threads", po::value<unsigned int>(&n_threads)->default_value(1), "Number of threads to use in computation")
		("incl-markers", po::value<vector<string> >()->composing(), "Comma-separated list of markers to include")
		("one-sided", po::bool_switch(&_onesided), "Generate pairwise models with one side given by the list of included markers or traits")
		;

	opts.add(regress_opts);

	return opts;
}

void Regression::printVarHeader(const string& var_name){
	out_f << var_name << "_Pval" << sep
		  << var_name << "_beta" << sep
		  << var_name << "_SE" << sep;
}

void Regression::printMarkerHeader(const string& var_name){
	if(encoding == Encoding::CODOMINANT){
		printVarHeader(var_name + "_Het");
		printVarHeader(var_name + "_Hom");
	}else{
		printVarHeader(var_name);
	}
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

	if(vm.count("excl-traits")){
		InputManager::parseInput(vm["excl-traits"].as<vector<string> >(), excl_traits);
	}

	if(vm.count("incl-markers")){
		InputManager::parseInput(vm["incl-markers"].as<vector<string> >(), incl_marker_name);
	}


	out_f.open(out_fn.c_str());
	if(!out_f){
		throw std::logic_error("Cannot open regression output file '" + out_fn + "'");
	}

	_threaded = n_threads > 0;
	if(model_files.size() == 0){
		if(_onesided && !pairwise){
			Logger::log_err("WARNING: --one-sided must be used with --pairwise; ignoring --one-sided directive");
			_onesided = false;
		}

		if(exclude_markers && !include_traits){
			Logger::log_err("ERROR: --exclude-markers requires --use-traits", true);
		}
	}else if (include_traits || _onesided || pairwise){
		Logger::log_err("WARNING: --use-traits, --one-sided, and --pairwise have no effect when specifying models");
		include_traits = _onesided = pairwise = false;
	}

	//OK, let's actually open these model files and read them!!
	vector<string>::const_iterator mf_itr = model_files.begin();
	while (mf_itr != model_files.end()) {

		ifstream input((*mf_itr).c_str());

		if (!input.is_open()) {
			Logger::log_err("ERROR: Error opening model file: " + (*mf_itr), true);
		}

		string line;
		while (getline(input, line)) {
			// string isn't empty or starts with "#"
			boost::algorithm::trim(line);
			if(line.size() != 0 && line[0] != '#'){
				_models.push_back(line);
			}
		}

		input.close();
		++mf_itr;
	}

}

void Regression::runRegression(const DataSet& ds){
	set<string> all_traits;
	if(_onesided){
		all_traits.insert(ds.beginTrait(), ds.endTrait());
		// remove excluded traits
		set<string> tmp_alltrait;
		std::set_difference(all_traits.begin(), all_traits.end(),
						    excl_traits.begin(), excl_traits.end(),
						    std::inserter(tmp_alltrait, tmp_alltrait.begin()));

		all_traits.clear();

		// remove covariates
		std::set_difference(tmp_alltrait.begin(), tmp_alltrait.end(),
				            covar_names.begin(), covar_names.end(),
				            std::inserter(all_traits, all_traits.begin()));

		// remove the outcome variable
		set<string>::iterator alltrait_out_itr = all_traits.find(outcome_name);
		if(alltrait_out_itr != all_traits.end()){
			all_traits.erase(alltrait_out_itr);
		}
	}

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

		// And finally, remove the outcome variable
		tmp_itr = incl_traits.find(outcome_name);
		if(tmp_itr != incl_traits.end()){
			incl_traits.erase(tmp_itr);
		}

	}else{
		incl_traits.clear();
	}

	// Get the set of markers to include, if needed
	set<const Marker*> marker_incl;
	set<string>::const_iterator mni = incl_marker_name.begin();
	while(mni != incl_marker_name.end()){
		const Marker* found_marker = ds.getMarker(*mni);
		if(found_marker == 0){
			Logger::log_err("WARNING: could not find marker: " + *mni + ", ignoring");
		}else{
			marker_incl.insert(found_marker);
		}
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

	if(outcome_name.size() > 0 && !ds.isTrait(outcome_name)){
		Utility::Logger::log_err("ERROR: Outcome '" + outcome_name + "' is not a recognized trait, aborting.", true);
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
		n_trait = m->traits.size();
		delete m;
	} else{
		// Calculate the number of SNPs / Env vars based on the options passed
		// in to the Regression object.
		n_snp = (!exclude_markers) * (1 + pairwise);
		n_trait = (incl_traits.size() > 0) * (1 + exclude_markers);
	}

	printHeader(n_snp, n_trait);

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
		++si;
	}

	this->initData(model_str, ds);

	ModelGenerator* mgp;

	if(_models.size() == 0){
		if(_onesided){
			mgp = new OneSidedModelGenerator(ds, marker_incl, incl_traits, all_traits, exclude_markers);
		} else if(marker_incl.size() > 0){
			mgp = new BasicModelGenerator<set<const Marker*>::const_iterator>(ds, marker_incl.begin(), marker_incl.end(), incl_traits, pairwise, exclude_markers);
		} else {
			mgp = new BasicModelGenerator<DataSet::const_marker_iterator>(ds, ds.beginMarker(), ds.endMarker(), incl_traits, pairwise, exclude_markers);
		}
	}else{
		mgp = new TargetedModelGenerator(ds, _models);
	}


	if (_threaded) {
		boost::thread_group all_threads;

		for (unsigned int i = 0; i < n_threads; i++) {
			boost::thread* t = new boost::thread(&Regression::start,this,boost::ref(*mgp), boost::cref(ds));
			all_threads.add_thread(t);

			//all_threads.create_thread(boost::bind(&Regression::start, boost::ref(this),boost::ref(mg), boost::cref(ds)));
		}
		all_threads.join_all();

	} else {
		// go here for debugging purposes, i.e. set --threads to 0
		start(*mgp, ds);
	}

	delete mgp;

	printResults();
}

void Regression::start(ModelGenerator& mg, const DataSet& ds){

	Model* nm;

	// synchronize
	_model_gen_mutex.lock();
	unsigned int lc=1;
	nm = mg();
	_model_gen_mutex.unlock();
	// end synchronize

	while( nm ){
		Result* r = run(nm, ds);

		// Run some univariate models
		if(show_uni && nm->markers.size() + nm->traits.size() > 1){
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

		// synchronize
		_result_mutex.lock();
		results.push_back(r);
		_result_mutex.unlock();
		// end synchronize

		delete nm;

		// give someone else a chance ...
		boost::this_thread::yield();

		// synchronize
		_model_gen_mutex.lock();
		++lc;
		nm = mg();
		_model_gen_mutex.unlock();
		// end synchronize
	}
}

void Regression::printHeader(unsigned int n_snp, unsigned int n_trait) {
	// Now, print the header
	for (unsigned int i = 0; i < n_snp; i++) {
		out_f << "Var" << i + 1 << "_ID" << sep << "Var" << i + 1 << "_Pos"
				<< sep << "Var" << i + 1 << "_MAF" << sep;
	}

	for (unsigned int i = 0; i < n_trait; i++) {
		out_f << "Var" << n_snp + i + 1 << "_ID" << sep;
	}

	out_f << "N_Missing" << sep;

	printExtraHeader();

	if (n_snp + n_trait > 1) {
		for (unsigned int i = 0; show_uni && i < n_snp + n_trait; i++) {
			string hdr = "Uni_Var" + boost::lexical_cast<string>(i + 1);
			if (i < n_snp) {
				printMarkerHeader(hdr);
			} else {
				printVarHeader(hdr);
			}
		}

		if (interactions) {
			for (unsigned int i = 0; i < n_snp + n_trait; i++) {
				string hdr = "Red_Var" + boost::lexical_cast<string>(i + 1);
				if (i < n_snp) {
					printMarkerHeader(hdr);
				} else {
					printVarHeader(hdr);
				}
			}

			for (unsigned int i = 0; i < n_snp + n_trait; i++) {
				string hdr = "Full_Var" + boost::lexical_cast<string>(i + 1);
				if (i < n_snp) {
					printMarkerHeader(hdr);
				} else {
					printVarHeader(hdr);
				}
			}

			for (unsigned int i = 0; i < n_snp + n_trait; i++) {
				for (unsigned int j = i + 1; j < n_snp + n_trait; j++) {
					string v1_hdr = "Full_Var" + boost::lexical_cast<string>(i
							+ 1);
					string v2_hdr = "_Var" + boost::lexical_cast<string>(j + 1);

					if (encoding == Encoding::CODOMINANT && i < n_snp) {

						if (j < n_snp) {
							printVarHeader(v1_hdr + "_Het" + v2_hdr + "_Het");
							printVarHeader(v1_hdr + "_Het" + v2_hdr + "_Hom");
							printVarHeader(v1_hdr + "_Hom" + v2_hdr + "_Het");
							printVarHeader(v1_hdr + "_Hom" + v2_hdr + "_Hom");

						} else {
							printVarHeader(v1_hdr + "_Het" + v2_hdr);
							printVarHeader(v1_hdr + "_Hom" + v2_hdr);
						}

					} else {
						printVarHeader(v1_hdr + v2_hdr);
					}
				}
			}

			out_f << "Red_Model_Pval" << sep;
		}
	}

	for (unsigned int i = 0; !interactions && i < n_snp + n_trait; i++) {
		string hdr = "Var" + boost::lexical_cast<string>(i + 1);
		if (i < n_snp) {
			printMarkerHeader(hdr);
		} else {
			printVarHeader(hdr);
		}
	}

	// If we are looking at single SNP models w/ categorical weight, print said weight!
	if (encoding == Encoding::WEIGHTED && n_snp == 1 && n_trait == 0) {
		out_f << "Categ_Weight" << sep;
	}

	out_f << "Overall_Pval";

	set<CorrectionModel>::const_iterator c_itr = corr_methods.begin();
	while (c_itr != corr_methods.end()) {
		out_f << sep << "Overall_Pval_adj_" << *c_itr;
		++c_itr;
	}

	out_f << std::endl;
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
	mod.categorical = true;
	Result* r = run(&mod, ds);

	// NOTE: this assigns to the map at the same step
	float toret = (categ_weight[m] = r->coeffs[0] / (r->coeffs[0] + r->coeffs[1]));

	delete r;

	return toret;

}

Regression::Result* Regression::run(const Model* m, const DataSet& ds) {

	unsigned int numLoci = m->markers.size();
	unsigned int numCovars = covar_names.size();
	unsigned int numTraits = m->traits.size();
	unsigned int n_vars = numLoci + numTraits;
	unsigned int n_interact = 0;

	// determine size of row for each sample in dataset
	unsigned int n_cols = 1 + numCovars + numLoci + numTraits ;
	if(m->categorical){
		n_cols += numLoci - numTraits;
	} else if (interactions){
		n_interact += (n_vars * (n_vars - 1))/2;
		n_cols += n_interact;
	}

	if(encoding == Encoding::CODOMINANT){
		// If we have codominant encoding, we need an extra column for every SNP
		n_cols += numLoci;

		// If we have interactions (bleh!), we need an extra column for every
		// SNP-Trait pair and 3 extra columns for every SNP-SNP pair!
		if(interactions){
			unsigned int toadd = numLoci * numTraits + 3 * numLoci * (numLoci - 1) / 2;
			n_interact += toadd;
			n_cols += toadd;
		}


	}

	// Allocate a huge amount of memory for the regression here
	double* regress_output = new double[_pheno.size()];

	double* regress_buf = new double[_pheno.size()* n_cols];

	boost::multi_array_ref<double, 2> regress_data(regress_buf, boost::extents[_pheno.size()][n_cols]);
	//typedef boost::multi_array<double, 3> array_type;
	//typedef array_type::index index;
	//  array_type A(boost::extents[3][4][2])

	double* row_data = new double[n_cols];
	unsigned char geno[numLoci];
	float maf_sum[numLoci];

	DataSet::const_sample_iterator si = ds.beginSample();
	DataSet::const_sample_iterator se = ds.endSample();

	unsigned int n_samples = 0;
	unsigned int n_missing = 0;

	for(unsigned int i=0; i<numLoci; i++){
		maf_sum[i] = 0;
	}

	double categ_weight[numLoci];
	if ((!m->categorical) && encoding == Encoding::WEIGHTED){
		for (unsigned int i = 0; i<numLoci; i++){
			categ_weight[i] = getCategoricalWeight(m->markers[i], ds);
		}
	}
	//DataSet::const_sample_iterator se = ds.endSample();

	while(si != se){
		unsigned int pos=0;
		// 1st col is always 1
		row_data[pos++] = 1;

		// Now, the covariates
		for(unsigned int i=0; i<_covars.size(); i++){
			row_data[pos++] = _covars[i][n_samples];
		}

		// Now, work through the markers
		for (unsigned int i = 0; i<numLoci; i++){
			geno[i] = (*si)->getAdditiveGeno(*(m->markers[i]));

			if(geno[i] == Sample::missing_allele){
				row_data[pos++] = numeric_limits<double>::quiet_NaN();
				if(m->categorical || encoding == Encoding::CODOMINANT){
					row_data[pos++] = numeric_limits<double>::quiet_NaN();
				}
			} else {
				if(m->categorical){
					row_data[pos++] = EncodingModel(Encoding::DOMINANT)(geno[i]);
					row_data[pos++] = EncodingModel(Encoding::RECESSIVE)(geno[i]);
				} else if (encoding == Encoding::CODOMINANT){
					row_data[pos++] = (geno[i] == 1);
					row_data[pos++] = (geno[i] == 2);
				}else if (encoding == Encoding::WEIGHTED){

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
		for(unsigned int i = 0; (!m->categorical) && i < numTraits; i++){
			row_data[pos++] = ds.getTrait(m->traits[i], *si);
		}

		// Now, the interaction terms
		for (unsigned int i = 0; interactions && (!m->categorical) && i	< (n_vars); i++) {
			for (unsigned int j = i + 1; j < (n_vars); j++) {

				if(encoding == Encoding::CODOMINANT){
					// we have to be a little careful here.
					if(i < numLoci){
						if(j < numLoci){
							row_data[pos++] = row_data[i*2 + numCovars + 1] * row_data[j*2 + numCovars + 1];
							row_data[pos++] = row_data[i*2 + numCovars + 1] * row_data[j*2 + numCovars + 2];
							row_data[pos++] = row_data[i*2 + numCovars + 2] * row_data[j*2 + numCovars + 1];
							row_data[pos++] = row_data[i*2 + numCovars + 2] * row_data[j*2 + numCovars + 2];
						} else {
							row_data[pos++] = row_data[i*2 + numCovars + 1] * row_data[j + numLoci + numCovars + 1];
							row_data[pos++] = row_data[i*2 + numCovars + 2] * row_data[j + numLoci + numCovars + 1];
						}
					} else {
						row_data[pos++] = row_data[i + numLoci + numCovars + 1] * row_data[j + numLoci + numCovars + 1];
					}

				}else{
					row_data[pos++] = row_data[i + numCovars + 1] * row_data[j + numCovars + 1];
				}
			}
		}

		// OK, check for missingness before adding it to the dataset
		bool ismissing = std::isnan(_pheno[n_samples]);
		for(unsigned int i=0; i<pos; i++){
			ismissing |= std::isnan(row_data[i]);
		}

		if(!ismissing){
			for(unsigned int i=0; i<numLoci; i++){
				maf_sum[i] += geno[i];
			}
			// do a fast memory copy into the data for regression
			std::memcpy(&regress_data[n_samples - n_missing][0], row_data, sizeof(double) * (pos));
			regress_output[n_samples - n_missing] = _pheno[n_samples];
		} else {
			++n_missing;
		}

		++n_samples;
		++si;
	}

	// The number of variables in the "reduced" model is:
	// # of covariates if no interactions
	// # of main effects (total columns - # interactions) o/w
	unsigned int red_vars = n_interact == 0 ? numCovars : n_cols - n_interact - 1;

	stringstream ss;

	for(unsigned int i=0; i<numLoci; i++){
		float maf = maf_sum[i] / (2*static_cast<float>(n_samples-n_missing));
		ss << m->markers[i]->getID() << sep
		   << m->markers[i]->getChromStr() << ":" << m->markers[i]->getLoc()
		   << sep << std::min(maf, 1-maf) << sep;
	}
	for(unsigned int i=0; i<numTraits; i++){
		ss << m->traits[i] << sep;
	}

	Result* r = 0;

	if(n_cols <= n_samples - n_missing){
		r = calculate(regress_output, &regress_data[0][0],
					n_cols, n_samples - n_missing, 0, red_vars);
	}else{
		Logger::log_err("WARNING: not enough samples in model: '" + ss.str() + "'!");
		r = new Result();
		r->p_val = 1;
		r->log_likelihood = r->r_squared = std::numeric_limits<float>::quiet_NaN();
	}

	// print the # missing from this model
	ss << n_missing << sep;

	r->prefix = ss.str();

	ss.clear();
	if(encoding == Encoding::WEIGHTED && m->markers.size() == 1 && m->traits.size() == 0 && !m->categorical){
		//ss << categ_weight[0] << sep;
		r->suffix += boost::lexical_cast<string>(categ_weight[0]) + sep;
	}

	delete[] regress_output;
	delete[] regress_buf;
	delete[] row_data;

	return r;
}

void Regression::addResult(Result* curr_result, const Result* null_result){
	if(null_result){
		unsigned int s = null_result->coeffs.size();

		for(int i=s-1; i >= 0; --i){
			curr_result->coeffs.push_front(null_result->coeffs[i]);
			curr_result->stderr.push_front(null_result->stderr[i]);
			curr_result->p_vals.push_front(null_result->p_vals[i]);
		}

		// If we have coefficients, it means that the null model is a
		// model with explanatory variables, and we need to add a suffix containing
		// the reduced model's p-value
		if(s > 0){
			stringstream ss;
			ss << null_result->p_val << sep;

			curr_result->suffix = ss.str();
		}
	}
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
			++c_itr;
		}

	}


	for(unsigned int i=0; i<results.size(); i++){
		// If we've gone over the threshold, stop printing!
		if(results[i]->p_val > cutoff_p){
			break;
		}
		// print a single line
		out_f << results[i]->prefix;

		printExtraResults(*(results[i]));

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
			++p_itr;
		}

		out_f << std::endl;
	}
}

Regression::Model* Regression::TargetedModelGenerator::next() {

	Model* m = 0;
	// models are given one per line here
	if(_mitr != _mend){
		m = Regression::parseModelStr(*_mitr, _ds);
		++_mitr;
	}

	return m;
}

Regression::Model* Regression::OneSidedModelGenerator::next() {

	Model* m = 0;

	// We want Env vars
	if (_traits) {
		// we want exhaustive EnvxEnv
		if (_nomarker) {
			if (_ti2 == _tall_set.end() && _titr != _t_set.end()) {
				_t_processed.insert(*(_titr++));
				if (_titr != _t_set.end()) {
					_ti2 = _tall_set.begin();
				}
			}

			// If ti2 == end ==> _titr == end
			if (_ti2 != _tall_set.end()) {
				if (*_titr != *_ti2 && _t_processed.find(*_ti2) == _t_processed.end()) {
					m = new Model();
					m->traits.push_back(*_titr);
					m->traits.push_back(*_ti2);
					++_ti2;
				} else {
					//If we're here, we have already seen this model!
					++_ti2;
					m = next();
				}
			}

			// We want exhaustive SNPxSNPxEnv models
		} else {
			// A little change from the BasicModelGenerator - here we
			// want to iterate over mi1, then mi2, then traits
			// (Basic iterated over traits, then mi1, then mi2)

			if (_titr == _t_set.end()) {
				// reset the traits and increment mi2
				if (_mi2 != _ds.endMarker() && ++_mi2 == _ds.endMarker()) {
					_titr = _t_set.begin();
					// OK, increment mi1 and add it to the "processed" list
					if (_mi1 != _m_set.end()) {
						_m_processed.insert(*(_mi1++));

						if (_mi1 != _m_set.end()) {
							_mi2 = _ds.beginMarker();
						} else {
							// If I'm here, we're done iterating!
							_titr = _t_set.end();
						}
					}
				}
			}

			if (_mi1 != _m_set.end()) {
				if (*_mi1 != *_mi2 && _m_processed.find(*_mi2) == _m_processed.end()) {
					m = new Model();
					m->markers.push_back(*_mi1);
					m->markers.push_back(*_mi2);
					m->traits.push_back(*_titr);
					++_mi2;
				} else {
					++_mi2;
					m = next();
				}

			}
		}
		// OK, we just want SNPxSNP models
	} else {

		if (_mi2 == _ds.endMarker() && _mi1 != _m_set.end()) {
			_m_processed.insert(*(_mi1++));
			if(_mi1 != _m_set.end()){
				_mi2 = _ds.beginMarker();
			}
		}

		if (_mi2 != _ds.beginMarker()) {
			if(*_mi2 != *_mi1 && _m_processed.find(*_mi2) == _m_processed.end()){
				m = new Model();
				m->markers.push_back(*_mi1);
				m->markers.push_back(*_mi2);
				++_mi2;
			} else {
				++_mi2;
				m = next();
			}
		}

	}

	// NOTE: m == 0 if no more models can be generated!
	// IN targeted mode, m != 0, but m->markers.size() == 0 && m->traits.size() == 0 in the case of a "bad" model
	return m;

}

}
}

