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
#include <cstdio>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include <boost/iostreams/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

//#define binary_iarchive text_iarchive
//#define binary_oarchive text_oarchive

#include <boost/serialization/export.hpp>


#include "config.h"
#ifdef HAVE_CXX_MPI
#include <mpi.h>
#endif

using std::string;
using std::vector;
using std::set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::min;
using std::max;
using std::ifstream;
using std::deque;
using std::multimap;
using std::pair;
using std::make_pair;

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Data::Sample;
using PLATO::Utility::InputManager;
using PLATO::Utility::Logger;

namespace po=boost::program_options;

BOOST_CLASS_EXPORT(PLATO::Analysis::Regression::ExtraData)

namespace PLATO{

namespace Analysis{

//BOOST_CLASS_EXPORT_GUID(Regression::ExtraData,"red")

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
	if(_lowmem){
		tmp_f.close();
	}

	if(mgp){
		delete mgp;
	}

	if(class_data){
		delete class_data;
	}
}

po::options_description& Regression::addOptions(po::options_description& opts){
	po::options_description model_opts("Model Generation Options");

	model_opts.add_options()
		("models", po::value<vector<string> >(&model_files)->composing(), "List of files containing models to generate")
		("exclude-markers", po::bool_switch(&exclude_markers),
			"Do not include markers in the generated models (use for EWAS)")
		("use-traits", po::bool_switch(&include_traits),
			"Include the traits not used as covariates in the generated models")
		("pairwise", po::bool_switch(&pairwise),
			"When auto-generating models, include two variables exhaustively")
		("incl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to include")
		("excl-traits", po::value<vector<string> >()->composing(), "Comma-separated list of traits to exclude")
		("incl-markers", po::value<vector<string> >()->composing(), "Comma-separated list of markers to include")
		("one-sided", po::bool_switch(&_onesided), "Generate pairwise models with one side given by the list of included markers or traits")
		;

	po::options_description regress_opts("Regression Options");

	regress_opts.add_options()
		("interactions", po::bool_switch(&interactions), "Include interactions in the models generated")
		("covariates", po::value<vector<string> >()->composing(), "A list of covariates to use in the model")
		("outcome", po::value<vector<string> >()->composing(), "Use a given covariate as the regression model outcome")
		("encoding", po::value<EncodingModel>(&encoding)->default_value("additive"), "Encoding model to use in the regression (additive, dominant, recessive, weighted, codominant)")
		("show-univariate", po::bool_switch(&show_uni), "Show univariate results in multivariate models")
		("phewas", po::bool_switch(&_phewas), "Perform a pheWAS (use all traits not included as covariates or specifically included)")
		("correction", po::value<vector<string> >()->composing(), ("p-value correction method(s) (" + Correction::listCorrectionMethods() + ")").c_str())
		("permutations", po::value<unsigned int>(&n_perms)->default_value(0), "Number of permutations to use in permutation testing (disabled by default - set to 0 to disable permutation)")
		("thresh", po::value<float>(&cutoff_p)->default_value(1.0f), "Threshold for printing resultant models")
		("output", po::value<string>(&out_fn)->default_value("output.txt"), "Name of the file to output results")
		("separator", po::value<string>(&sep)->default_value("\t", "<TAB>"), "Separator to use when outputting results file")
		("threads", po::value<unsigned int>(&n_threads)->default_value(1), "Number of threads to use in computation")
		("lowmem", po::bool_switch(&_lowmem), "Reduce the memory footprint (at a potential performance penalty)")
		;

	opts.add(model_opts).add(regress_opts);

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

	if(vm.count("outcome")){
		InputManager::parseInput(vm["outcome"].as<vector<string> >(), outcome_names);
	}

	if(_phewas){
		if(vm.count("outcome")){
			Logger::log_err("WARNING: --outcome and --phewas are incompatible, ignoring --phewas");
			_phewas = false;
		} else if(include_traits && !vm.count("incl-traits")){
			Logger::log_err("WARNING: --phewas used with --use-traits and without --incl-traits, ignoring --phewas");
			_phewas = false;
		}


	}

	out_f.open(out_fn.c_str());
	if(!out_f){
		throw std::logic_error("Cannot open regression output file '" + out_fn + "'");
	}

#ifdef HAVE_CXX_MPI
	int n_procs = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	if(n_procs > 1){
		_use_mpi = true;
	}
#endif

	if(_use_mpi){
		if(n_threads > 1){
			Logger::log_err("WARNING: --threads used with MPI, setting number of threads to 0");
		}
		n_threads = 0;
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

	if(_lowmem){
		tmp_f.open(boost::iostreams::file_descriptor(fileno(std::tmpfile())),
				std::ios_base::binary | std::ios_base::in | std::ios_base::out);
	}

}

void Regression::runRegression(const DataSet& ds){

	if(_phewas){
		// I want the outcome_names to consist of all (enabled) traits, except for those
		// specifically excluded, included (in the case of use_traits), or used as covariates

		// outcome_names should be empty at this point, but let's make sure
		outcome_names.clear();

		outcome_names.insert(ds.beginTrait(), ds.endTrait());

		set<string> outcome_tmp;
		std::set_difference(outcome_names.begin(), outcome_names.end(),
				            excl_traits.begin(), excl_traits.end(),
				            std::inserter(outcome_tmp, outcome_tmp.begin()));

		outcome_names.clear();

		std::set_difference(outcome_tmp.begin(), outcome_tmp.end(),
						    covar_names.begin(), covar_names.end(),
						    std::inserter(outcome_names, outcome_names.begin()));

		if(include_traits){
			outcome_tmp.clear();
			std::set_difference(outcome_names.begin(), outcome_names.end(),
					            incl_traits.begin(), incl_traits.end(),
					            std::inserter(outcome_tmp, outcome_tmp.begin()));

			outcome_names.clear();
			outcome_names.insert(outcome_tmp.begin(), outcome_tmp.end());
		}

	}

	// Make sure all the outcomes are actual traits
	set<string>::iterator out_itr = outcome_names.begin();
	while(out_itr != outcome_names.end()){
		if(!ds.isTrait(*out_itr)){
			Logger::log_err("WARNING: '" + *out_itr + "' is not a recognized covariate, ignoring.");
			outcome_names.erase(out_itr++);
		} else{
			++out_itr;
		}
	}


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

		tmp_alltrait.clear();
		// remove the outcome variable(s)
		std::set_difference(all_traits.begin(), all_traits.end(),
							outcome_names.begin(), outcome_names.end(),
							std::inserter(tmp_alltrait, tmp_alltrait.begin()));

		all_traits.clear();
		all_traits.insert(tmp_alltrait.begin(), tmp_alltrait.end());
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

		// now, take out the outcomes
		tmp_trait.clear();
		std::set_difference(incl_traits.begin(), incl_traits.end(),
						    outcome_names.begin(), outcome_names.end(),
						    std::inserter(tmp_trait, tmp_trait.begin()));

		incl_traits.clear();
		incl_traits.insert(tmp_trait.begin(), tmp_trait.end());

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

	unsigned int n_snp = 0;
	unsigned int n_trait = 0;

	if(_models.size() > 0){
		deque<string>::const_iterator mi = _models.begin();
		Model* m = 0;
		while(mi != _models.end() && (m ==0 || (m->traits.size() == 0 && m->markers.size() ==0))){
			m = parseModelStr(*mi, ds);
			++mi;
		}

		n_snp = m->markers.size();
		n_trait = m->traits.size();

		delete m;

		if(n_snp ==0 && n_trait == 0){
			Logger::log_err("ERROR: Could not find a model to run", true);
		}

	} else{
		// Calculate the number of SNPs / Env vars based on the options passed
		// in to the Regression object.
		n_snp = (!exclude_markers) * (1 + pairwise);
		n_trait = (incl_traits.size() > 0) * (1 + pairwise * exclude_markers);
	}

	// Add in  the extra dfs, along with their column IDs
	if(encoding == Encoding::WEIGHTED){
		for(unsigned int i=0; i<n_snp; i++){
			_extra_df_map[covar_names.size() + 1 + i] = 1;
		}
	}

	printHeader(n_snp, n_trait);

	DataSet::const_sample_iterator si = ds.beginSample();

	// set up the outcome vector
	if(outcome_names.size() == 0){
		while (si != ds.endSample()) {
			// This will be NaN for missing, 0/1 for case/control
			// or a quantitative value as appropriate
			_pheno.push_back((*si)->getPheno());
			++si;
		}
		outcome_names.insert("");
		// put it back!!
		si = ds.beginSample();

	}

	// Now, set up the outcome variable and the covariates
	while(si != ds.endSample()){
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

	// set up the model generator
	if (_models.size() == 0) {
		if (_onesided) {
			mgp = new OneSidedModelGenerator(ds, marker_incl,
					incl_traits, all_traits, exclude_markers);
		} else if (marker_incl.size() > 0) {
			mgp = new BasicModelGenerator<
					set<const Marker*>::const_iterator> (ds,
					marker_incl.begin(), marker_incl.end(),
					incl_traits, pairwise, exclude_markers);
		} else {
			mgp = new BasicModelGenerator<
					DataSet::const_marker_iterator> (ds,
					ds.beginMarker(), ds.endMarker(), incl_traits,
					pairwise, exclude_markers);
		}
	} else {
		mgp = new TargetedModelGenerator(ds, _models);
	}

	output_itr = outcome_names.begin();

	ds_ptr = &ds;

	if(*output_itr != "" || !initData()){
		if(*output_itr == ""){
			++output_itr;
		}
		while(!resetPheno(*output_itr) && ++output_itr != outcome_names.end());
	}

	if(_use_mpi){
		processMPI();
	} else {
		while(output_itr != outcome_names.end()){
			// set up the phenotype (if we haven't already!)
			si = ds.beginSample();

			mgp->reset();

			if (_threaded) {
				boost::thread_group all_threads;

				for (unsigned int i = 0; i < n_threads; i++) {
					boost::thread* t = new boost::thread(&Regression::start,
							this, boost::cref(ds),	boost::cref(*output_itr));
					all_threads.add_thread(t);

					//all_threads.create_thread(boost::bind(&Regression::start, boost::ref(this),boost::ref(mg), boost::cref(ds)));
				}
				all_threads.join_all();

			} else {
				// go here for debugging purposes, i.e. set --threads to 0
				start(ds, *output_itr);
			}

			while(++output_itr != outcome_names.end() && !resetPheno(*output_itr));

		}
	}

	printResults();
}

void Regression::start(const DataSet& ds, const string& outcome){

	Model* nm;

	// synchronize
	_model_gen_mutex.lock();
	unsigned int lc=1;
	nm = mgp->next();
	_model_gen_mutex.unlock();
	// end synchronize

	while( nm ){
		Result* r = run(*nm);

		// run will return 0 if something went wrong (i.e. did not run anything)
		if (r) {

			// Run some univariate models
			if (show_uni && nm->markers.size() + nm->traits.size() > 1) {
				addUnivariate(*r, *nm);
			}
			// put the outcome right at the beginning of the prefix!
			r->prefix = outcome + sep + r->prefix;

			addResult(r);
		}

		delete nm;

		// give someone else a chance ...
		boost::this_thread::yield();

		// synchronize
		_model_gen_mutex.lock();
		++lc;
		nm = mgp->next();
		_model_gen_mutex.unlock();
		// end synchronize
	}
}

void Regression::addResult(Result* r){
	if(r){
		_result_mutex.lock();
		if (_lowmem) {
			result_pvals.push_back(r->p_val);
			printResultLine(*r, tmp_f);
			tmp_f << "\n";
			delete r;
		} else {
			results.push_back(r);
		}
		_result_mutex.unlock();
	}

}

bool Regression::resetPheno(const string& pheno){

	_pheno.clear();
	DataSet::const_sample_iterator si = ds_ptr->beginSample();
	while(si != ds_ptr->endSample()){
		_pheno.push_back(ds_ptr->getTrait(pheno,*si));
		++si;
	}

	// clear the data structures that depend on the phenotype
	_marker_uni_result.clear();
	_trait_uni_result.clear();
	categ_weight.clear();

	// if we have a phenotype that is not appropriate, move on to the
	// next one

	return initData();
}

void Regression::printHeader(unsigned int n_snp, unsigned int n_trait) {
	// Now, print the header
	out_f << "Outcome" << sep;

	for (unsigned int i = 0; i < n_snp; i++) {
		out_f << "Var" << i + 1 << "_ID" << sep << "Var" << i + 1 << "_Pos"
				<< sep << "Var" << i + 1 << "_MAF" << sep;
	}

	for (unsigned int i = 0; i < n_trait; i++) {
		out_f << "Var" << n_snp + i + 1 << "_ID" << sep;
	}

	out_f << "Num_Missing" << sep;

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

			// These will come in the suffix of the result!
			out_f << "Red_Model_Pval" << sep << "Full_Model_Pval" << sep;
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
		const Marker* m = ds.getMarker(model_elements[i]);
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

void Regression::addUnivariate(Result& r, const Model& m){
	for (unsigned int i = 0; i < m.markers.size(); i++) {
		_univar_mmutex.lock();
		map<const Marker*, Result*>::const_iterator m_itr =
				_marker_uni_result.find(m.markers[i]);
		// If this is the case, we have not seen this marker before
		if (m_itr == _marker_uni_result.end()) {
			Model uni_m;
			uni_m.markers.push_back(m.markers[i]);
			m_itr = _marker_uni_result.insert(
					_marker_uni_result.begin(),
					std::make_pair(m.markers[i], run(uni_m)));
		}
		_univar_mmutex.unlock();

		r.unimodel.push_back((*m_itr).second);
	}

	for (unsigned int i = 0; i < m.traits.size(); i++) {
		_univar_tmutex.lock();
		map<string, Result*>::const_iterator t_itr =
				_trait_uni_result.find(m.traits[i]);
		// If this is the case, we have not seen this marker before
		if (t_itr == _trait_uni_result.end()) {
			Model m;
			m.traits.push_back(m.traits[i]);
			t_itr = _trait_uni_result.insert(
					_trait_uni_result.begin(),
					std::make_pair(m.traits[i], run(m)));
		}
		_univar_tmutex.unlock();

		r.unimodel.push_back((*t_itr).second);
	}
}

float Regression::getCategoricalWeight(const Marker* m){

	float toret = 0.5;

	_categ_mutex.lock();
	map<const Marker*, float>::const_iterator w_itr = categ_weight.find(m);
	if(w_itr != categ_weight.end()){
		toret = (*w_itr).second;
	} else {

		// OK, at this point, we know we don't have the categorical value

		Model mod;
		mod.markers.push_back(m);
		mod.categorical = true;
		Result* r = run(mod);

		// NOTE: this assigns to the map at the same step
		// Make sure to check for existence of the result!
		// If nonexistent (you probably have bigger problems), default to 1/2 (additive encoding)
		toret = (r ? r->coeffs[0] / (r->coeffs[0] + r->coeffs[1]) : 0.5);
		toret = categ_weight[m] = (std::isfinite(toret) ? toret : 0.5);

		if(r){
			delete r;
		}
	}
	_categ_mutex.unlock();

	return toret;

}

Regression::calc_matrix* Regression::getCalcMatrix(const Model& m){
	calc_matrix* ret_val = new calc_matrix();

	unsigned int numLoci = m.markers.size();
	unsigned int numCovars = covar_names.size();
	unsigned int numTraits = m.traits.size();
	unsigned int n_vars = numLoci + numTraits;
	unsigned int n_interact = 0;

	// determine size of row for each sample in dataset
	ret_val->n_cols = 1 + numCovars + numLoci + numTraits ;
	if(m.categorical){
		ret_val->n_cols += numLoci - numTraits;
	} else if (interactions){
		n_interact += (n_vars * (n_vars - 1))/2;
		ret_val->n_cols += n_interact;
	}

	if(encoding == Encoding::CODOMINANT){
		// If we have codominant encoding, we need an extra column for every SNP
		ret_val->n_cols += numLoci;

		// If we have interactions (bleh!), we need an extra column for every
		// SNP-Trait pair and 3 extra columns for every SNP-SNP pair!
		if(interactions){
			unsigned int toadd = numLoci * numTraits + 3 * numLoci * (numLoci - 1) / 2;
			n_interact += toadd;
			ret_val->n_cols += toadd;
		}
	}

	// Allocate a huge amount of memory for the regression here
	ret_val->outcome = new double[_pheno.size()];
	ret_val->data = new double[_pheno.size()* ret_val->n_cols];

	boost::multi_array_ref<double, 2> regress_data(ret_val->data, boost::extents[_pheno.size()][ret_val->n_cols]);
	//typedef boost::multi_array<double, 3> array_type;
	//typedef array_type::index index;
	//  array_type A(boost::extents[3][4][2])

	double row_data[ret_val->n_cols];
	unsigned char geno[numLoci];
	float maf_sum[numLoci];
	double categ_weight[numLoci];

	DataSet::const_sample_iterator si = ds_ptr->beginSample();
	DataSet::const_sample_iterator se = ds_ptr->endSample();

	unsigned int n_samples = 0;
	unsigned int n_missing = 0;

	for(unsigned int i=0; i<numLoci; i++){
		maf_sum[i] = 0;
	}

	if ((!m.categorical) && encoding == Encoding::WEIGHTED){
		for (unsigned int i = 0; i<numLoci; i++){
			categ_weight[i] = getCategoricalWeight(m.markers[i]);
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
			geno[i] = (*si)->getAdditiveGeno(*(m.markers[i]));

			if(geno[i] == Sample::missing_allele){
				row_data[pos++] = numeric_limits<double>::quiet_NaN();
				if(m.categorical || encoding == Encoding::CODOMINANT){
					row_data[pos++] = numeric_limits<double>::quiet_NaN();
				}
			} else {
				if(m.categorical){
					row_data[pos++] = EncodingModel(Encoding::DOMINANT)(geno[i]);
					row_data[pos++] = EncodingModel(Encoding::RECESSIVE)(geno[i]);
				} else if (encoding == Encoding::CODOMINANT){
					row_data[pos++] = (geno[i] == 1);
					row_data[pos++] = (geno[i] == 2);
				}else if (encoding == Encoding::WEIGHTED){
					// works if w in [0,1]
					row_data[pos++] = categ_weight[i] * (geno[i] == 1) + (geno[i] == 2);
				} else {
					row_data[pos++] = encoding(geno[i]);
				}
			}
		}

		// Now, work through the traits
		for(unsigned int i = 0; (!m.categorical) && i < numTraits; i++){
			row_data[pos++] = ds_ptr->getTrait(m.traits[i], *si);
		}

		// Now, the interaction terms
		for (unsigned int i = 0; interactions && (!m.categorical) && i	< (n_vars); i++) {
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
			ret_val->outcome[n_samples - n_missing] = _pheno[n_samples];
		} else {
			++n_missing;
		}

		++n_samples;
		++si;
	}

	// The number of variables in the "reduced" model is:
	// # of covariates if no interactions
	// # of main effects (total columns - # interactions) o/w
	ret_val->red_vars = n_interact == 0 ? numCovars : ret_val->n_cols - n_interact - 1;
	ret_val->n_sampl = n_samples - n_missing;

	if(ret_val->n_sampl < ret_val->n_cols){
		// If this is the case, we CANNOT run this!
		delete ret_val;
		ret_val = 0;
	} else {
		// Print some summary data helpful for the result that we get when
		// calculating the matrix to send.
		stringstream ss;
		for(unsigned int i=0; i<numLoci; i++){
			float maf = maf_sum[i] / (2*static_cast<float>(n_samples-n_missing));
			ss << m.markers[i]->getID() << sep
			   << m.markers[i]->getChromStr() << ":" << m.markers[i]->getLoc()
			   << sep << m.markers[i]->getAltAllele() << ":" << maf << sep;
		}
		for(unsigned int i=0; i<numTraits; i++){
			ss << m.traits[i] << sep;
		}
		ss << n_missing << sep;
		ret_val->prefix = ss.str();
		ss.clear();
	}

	return ret_val;
}

Regression::Result* Regression::run(const Model& m) {

	calc_matrix* all_data = getCalcMatrix(m);

	Result* r = 0;

	calc_fn& calculate = (getCalcFn());

	if(all_data){
		r = calculate(all_data->outcome, all_data->data, all_data->n_cols,
				all_data->n_sampl, 0, all_data->red_vars, getExtraData());
	}else{
		Logger::log_err("WARNING: not enough samples in model: '" + m.getID() + "'!");
		r = 0;
	}

	if(r){
		// print the # missing from this model
		r->prefix = all_data->prefix;

		if(encoding == Encoding::WEIGHTED && m.markers.size() == 1 && m.traits.size() == 0 && !m.categorical){
			//ss << categ_weight[0] << sep;
			r->suffix += boost::lexical_cast<string>(categ_weight[0]) + sep;
		}
	}

	if(all_data){
		delete all_data;
	}

	return r;
}

void Regression::printResults(){

	vector<std::iostream::pos_type> file_pos;
	vector<size_t> idx_pos;
	size_t n_results = _lowmem ? result_pvals.size() : results.size();
	string tmpf_line;

	if (_lowmem) {
		idx_pos.reserve(n_results);

		for(size_t i=0; i<n_results; i++){
			idx_pos.push_back(i);
		}

		std::sort(idx_pos.begin(), idx_pos.end(), pval_sorter(result_pvals));

		file_pos.reserve(n_results+1);
		// get the positions of the beginning of every line
		tmp_f.flush();
		tmp_f.seekg(0, std::ios_base::beg);
		file_pos.push_back(tmp_f.tellg());
		while( getline(tmp_f, tmpf_line) ){
			file_pos.push_back(tmp_f.tellg());
		}
		tmp_f.seekg(0, std::ios_base::beg);
		tmp_f.clear();
		//int res = boost::iostreams::seek(tmp_f, 0, std::ios_base::beg);
	} else {
		// Perhaps sort the deque of results based on overall p-values
		std::sort(results.begin(), results.end(), result_sorter());
	}

	// Now, let's do some multiple test correction!
	map<CorrectionModel, vector<float> > pval_corr;

	// don't bother with the work if there's no correction desired!
	if (corr_methods.size()) {

		vector<float> pv_in;
		if (_lowmem) {
			pv_in.reserve(n_results);
			for(size_t i=0; i<n_results; i++){
				pv_in.push_back(result_pvals[idx_pos[i]]);
			}
		} else {
			pv_in.reserve(n_results);
			for (size_t i = 0; i < n_results; i++) {
				pv_in.push_back(results[i]->p_val);
			}
		}

		set<CorrectionModel>::const_iterator c_itr = corr_methods.begin();
		while (c_itr != corr_methods.end()) {
			Correction::getCorrectionMethod(*c_itr)->correct(pv_in,	pval_corr[(*c_itr)]);
			++c_itr;
		}
	}

	for (size_t i = 0; i < n_results; i++) {
		// If we've gone over the threshold, stop printing!
		if ((_lowmem ? result_pvals[idx_pos[i]] : results[i]->p_val)  > cutoff_p) {
			break;
		}

		if(_lowmem){
			tmp_f.seekg(file_pos[idx_pos[i]], std::ios_base::beg);
			tmpf_line = "";
			std::getline(tmp_f, tmpf_line);
			out_f << tmpf_line;
			out_f << result_pvals[idx_pos[i]];
			tmp_f.clear();
		} else {
			printResultLine(*(results[i]), out_f);
			out_f << results[i]->p_val;
		}

		// print some creected p-values!
		// NOTE: we assume that "set" and "map" use the same ordering!!
		map<CorrectionModel, vector<float> >::const_iterator p_itr =
				pval_corr.begin();
		while (p_itr != pval_corr.end()) {
			out_f << sep << (*p_itr).second[i];
			++p_itr;
		}

		out_f << std::endl;
	}
}

Regression::Model* Regression::TargetedModelGenerator::next() {

	Model* m = 0;
	// models are given one per line here
	while(_mitr != _mend && (m == 0 || (m->markers.size() == 0 && m->traits.size() == 0))){
		m = Regression::parseModelStr(*_mitr, _ds);
		if(m->markers.size() == 0 && m->traits.size() == 0){
			Logger::log_err("WARNING: Model '" + *_mitr + "' has elements not in the dataset, ignoring");
			delete m;
			m = 0;
		}
		++_mitr;
	}

	return m;
}

void Regression::TargetedModelGenerator::reset() {
	_mitr = _mbegin;
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

void Regression::OneSidedModelGenerator::reset() {
	 //_m_set(m_set),
	//			_t_set(trait), _tall_set(all_traits),
	_mi1 = _m_set.begin();
	_mi2 = _ds.beginMarker();
	_titr = _t_set.begin();
	_ti2 = _tall_set.begin();
}

const Regression::ExtraData* Regression::getExtraData() const{
	if(!class_data){
		class_data = new ExtraData();
		class_data->encoding = encoding;
		class_data->base_covars = covar_names.size();
		class_data->sep = sep;
		class_data->n_extra_df_col = 0;
		for(map<unsigned int, unsigned int>::const_iterator dfi = _extra_df_map.begin();
				dfi != _extra_df_map.end(); dfi++){
			class_data->n_extra_df_col += (*dfi).second;
		}
		unsigned int* colptr = class_data->extra_df_col = new unsigned int[class_data->n_extra_df_col];
		for(map<unsigned int, unsigned int>::const_iterator dfi = _extra_df_map.begin();
				dfi != _extra_df_map.end(); dfi++){
			for(unsigned int i=0; i<(*dfi).second; i++){
				*(colptr++) = (*dfi).second;
			}
		}
	}

	return class_data;
}

unsigned int Regression::findDF(const gsl_matrix* P,
		unsigned int reduced_vars,
		unsigned int n_dropped,
		unsigned int* extra_cols,
		unsigned int n_extra_cols) {

	unsigned int n_cols = P->size1;
	unsigned int edf = n_cols - reduced_vars;

	// Now, we need to actually find the column indices that were kept, but only
	// if we dropped any
	if (n_dropped > 0) {
		gsl_vector* df_check = gsl_vector_calloc(n_cols);
		gsl_vector* df_check_t = gsl_vector_calloc(n_cols);
		for (unsigned int i = reduced_vars; i < n_cols; i++) {
			gsl_vector_set(df_check, i, 1);
		}

		// Now, permute the df_check
		gsl_blas_dgemv(CblasNoTrans, 1.0, P, df_check, 0.0, df_check_t);

		// unset the last # dropped
		for (unsigned int i = 1; i <= n_dropped; i++) {
			gsl_vector_set(df_check_t, n_cols - i, 0);
		}

		// unpermute
		gsl_blas_dgemv(CblasTrans, 1.0, P, df_check_t, 0.0, df_check);

		edf = gsl_blas_dasum(df_check);

		// Now, iterate over all of the elements in the extra_df_map:
		for (unsigned int i=0; i<n_extra_cols && extra_cols[i] < n_cols; i++){
			edf += gsl_vector_get(df_check, extra_cols[i]);
		}
		gsl_vector_free(df_check);
		gsl_vector_free(df_check_t);
	} else if(n_extra_cols > 0) {
		for (unsigned int i=0; i<n_extra_cols && extra_cols[i] < n_cols; i++){
			edf += (extra_cols[i] >= reduced_vars);
		}
	}

	return edf;
}

void Regression::printResult(const Result& r, std::ostream& of){
	if(r.submodel){
		printResult(*(r.submodel), of);
	}

	for(unsigned int j=0;j<r.n_vars; j++){
		of << r.p_vals[j] << sep
			  << r.coeffs[j] << sep
			  << r.stderr[j] << sep;
	}
}

void Regression::printResultLine(const Result& r, std::ostream& of){
	of << r.prefix;
	of << printExtraResults(r);

	for(unsigned int i=0; i<r.unimodel.size(); i++){
		if(r.unimodel[i]){
			printResult(*r.unimodel[i], of);
		} else{
			// If we're here, we couldn't run the univariate model
			// (likely something VERY bad has happened)
			of << numeric_limits<float>::quiet_NaN() << sep
			   << numeric_limits<float>::quiet_NaN() << sep
			   << numeric_limits<float>::quiet_NaN() << sep;
		}
	}

	printResult(r, of);
	if(r.submodel && r.submodel->n_vars > 0){\
		of << r.submodel->p_val << sep;
	}
	of << r.suffix;
}

pair<unsigned int, const char*> Regression::generateMsg(const Model& m){
	pair<unsigned int, const char*> retval(0,0);

	calc_matrix* data = getCalcMatrix(m);
	data->prefix = *output_itr + sep + data->prefix;

	// do the voodoo that you do so well here
	if(data){
		++msg_id;
		stringstream ss;
		boost::archive::binary_oarchive oa(ss);
		mpi_query q;
		q.msg_id = msg_id;
		q.calc_data = data;
		q.class_data = getExtraData();

		ss.flush();

		// save this to the archive
		oa << q;

		retval.first = ss.tellp();
		char* output = new char[retval.first];
		ss.read(output,retval.first);
		retval.second = output;

		work_map[msg_id] = &m;

		delete data;
	} else {
		// If we can't run this model, go back to square 1 (which will come back here, incidentially)
		retval = nextQuery();
	}

	return retval;
}

// MPI stuff here!
pair<unsigned int, const char*> Regression::nextQuery(){
	pair<unsigned int, const char*> retval(0,0);

	const Model * m;
	// If we have models queued up, please send them
	if(!model_queue.empty()){
		m = model_queue.front();
		model_queue.pop_front();
		retval = generateMsg(*m);
	// otherwise, generate a new model
	} else if( (m = mgp->next()) != 0) {
		// check for pre-locks (caused by weighted encoding)

		// NOTE: n_val_ptr will be 0 if no pre-locks were needed!
		int * n_val_ptr = 0;
		if(encoding == Encoding::WEIGHTED){

			for(unsigned int i=0; i<m->markers.size(); i++){
				if(categ_weight.find(m->markers[i]) == categ_weight.end()){
					Model *pm = new Model;
					pm->markers.push_back(m->markers[i]);
					pm->categorical = true;

					if(n_val_ptr == 0){
						n_val_ptr = new int(0);
						pre_lock_map.insert(make_pair(m, n_val_ptr));
					}

					++(*n_val_ptr);
					work_lock_map.insert(make_pair(pm->getID(), n_val_ptr));

					// only add if we aren't currently working on it!
					if(work_lock_map.count(pm->getID()) == 1){
						model_queue.push_back(pm);
					} else {
						delete pm;
					}

				}
			}
		}

		// if there ARE pre-locks, I will have added them to the queue
		if(n_val_ptr == 0){
			retval = generateMsg(*m);
		} else if(!model_queue.empty()){
			m = model_queue.front();
			model_queue.pop_front();
			retval = generateMsg(*m);
		} else  {
			retval = nextQuery();
		}
		// OK, now we can generate our message to send
		//retval = generateMsg(*m);
	} else if (pre_lock_map.size() != 0){
		// if I'm done with models, but I'm waiting on data, send that
		retval = pair<unsigned int, const char*>(pre_lock_map.size(), 0);
	} else if(show_uni && !work_map.empty()){
		// If I'm here, I might generate some post-locking models, so
		// send a "please wait" message
		retval = pair<unsigned int, const char*>(work_map.size(), 0);
	} else if(output_itr != outcome_names.end() && ++output_itr != outcome_names.end()){

		// Note the empty loop here - we will check and make sure that the
		// phenotype is good and continue until we either run out or find
		// a good one
		while(!resetPheno(*output_itr) && ++output_itr != outcome_names.end());

		// If I'm here, I have another phenotype to run!
		if(output_itr != outcome_names.end()){
			mgp->reset();
			// I don't want to recreate the logic above by creating a new model,
			// so just send a "we have more" signal and prepare to generate a new model
			retval = nextQuery();
		}
	}
	// empty else statement means that there is really NOTHING else to run!

	return retval;
}

void Regression::processResponse(unsigned int bufsz, const char* in_buf){
	// Turn this into a msg_id + Result

	stringstream ss_resp;
	ss_resp.write(in_buf, bufsz);
	ss_resp.seekg(0);
	boost::archive::binary_iarchive ia_resp(ss_resp);
	mpi_response resp;
	ia_resp >> resp;
	unsigned int msg = resp.msg_id;
	Result* r = resp.result;

	// Get the Model associated with the message
	const Model* m = 0;
	map<unsigned int, const Model*>::iterator witr = work_map.find(msg);
	if(witr != work_map.end()){
		m = witr->second;
		work_map.erase(witr);
	}

	// if this is found in the lock map, it's not a result!
	pair<multimap<string, int*>::iterator, multimap<string, int*>::iterator> lm_range = work_lock_map.equal_range(m->getID());
	bool is_result = (lm_range.first == lm_range.second);
	// decrement any locks we found!
	while(lm_range.first != lm_range.second){
		--(*((lm_range.first)->second));
		work_lock_map.erase(lm_range.first++);
	}

	// if this was a categorical test, add it as appropriate:
	if(m->categorical && m->markers.size() == 1){
		float wt = (r ? r->coeffs[0] / (r->coeffs[0] + r->coeffs[1]) : 0.5);
		categ_weight[m->markers[0]] = (std::isfinite(wt) ? wt : 0.5);
	} else if (!is_result &&  m->markers.size() +  m->traits.size() == 1){
		// if this is true, then the result should be added to the univariate
		// results
		if(m->markers.size() == 1){
			_marker_uni_result[m->markers[0]] = r;
		} else {
			_trait_uni_result[m->traits[0]] = r;
		}
	}

	// now that locks are decremented, see if anything needs to be placed
	// into the model queue.
	for(map<const Model*, int*>::iterator it = pre_lock_map.begin(); it != pre_lock_map.end(); ){
		if(*(it->second) == 0){
			delete it->second;
			model_queue.push_back(it->first);
			pre_lock_map.erase(it++);
		} else {
			++it;
		}
	}

	// Also, see if anything needs to be placed in the result list
	for(map<Result*, int*>::iterator it = post_lock_map.begin(); it != post_lock_map.end(); ){
		if(*(it->second) == 0){
			delete it->second;
			if(show_uni){
				map<Result*, const Model*>::iterator plm_itr = post_lock_models.find(it->first);
				if(plm_itr != post_lock_models.end()){
					//this really should always be the case!!  If not, we had
					// something very bad happen
					addUnivariate(*(it->first), *(plm_itr->second));
					// make sure to delete the model before erasing it from
					// the post_lock_model map!
					delete plm_itr->second;
					post_lock_models.erase(plm_itr);
				} else {
					Logger::log_err("WARNING: Unexpected error in displaying univariate models using MPI.");
				}
			}
			addResult(it->first);
			post_lock_map.erase(it++);
		} else {
			++it;
		}
	}

	vector<Model*> pl_needed;
	// check for post-locks needed, and add them to the queue
	if(show_uni && m->markers.size() +  m->traits.size() > 1){
		Model *uni_m;
		for (unsigned int i = 0; i < m->markers.size(); i++) {
			map<const Marker*, Result*>::const_iterator m_itr =
					_marker_uni_result.find(m->markers[i]);
			// If this is the case, we have not seen this marker before
			if (m_itr == _marker_uni_result.end()) {
				uni_m = new Model();
				uni_m->markers.push_back(m->markers[i]);
				pl_needed.push_back(uni_m);
			}
		}

		for (unsigned int i = 0; i < m->traits.size(); i++) {
			map<string, Result*>::const_iterator t_itr =
					_trait_uni_result.find(m->traits[i]);
			// If this is the case, we have not seen this marker before
			if (t_itr == _trait_uni_result.end()) {
				uni_m = new Model();
				uni_m->traits.push_back(m->traits[i]);
				pl_needed.push_back(uni_m);
			}
		}

		// in this case, it's safe to add the univariate results, b/c we
		// have already calculated them!
		if(pl_needed.size() == 0){
			addUnivariate(*r, *m);
		}
	}

	int* val_ptr = 0;
	for(unsigned int i=0; i<pl_needed.size(); i++){
		if(!val_ptr){
			val_ptr = new int(pl_needed.size());
			post_lock_map.insert(make_pair(r, val_ptr));
			post_lock_models.insert(make_pair(r,m));
		}
		// I know that the post-lock models can't need a pre-lock, phew!
		Model* plmod = pl_needed[i];
		work_lock_map.insert(make_pair(plmod->getID(), val_ptr));
		// Only add this to the model once, please!
		if(work_lock_map.count(plmod->getID()) == 1){
			model_queue.push_back(plmod);
		} else {
			delete plmod;
		}
		// if I needed post-locks, I'm not ready to push the result!
		is_result = false;
	}


	if(is_result){
		//is a result and I don't need any post-locks, add the result to the list
		addResult(r);
	}

	// now we're done with this model, so delete it, but only if we don't need
	// it later!
	if(pl_needed.size() == 0){
		delete m;
	}
}

pair<unsigned int, const char*> Regression::calculate_MPI(unsigned int bufsz, const char* in_buf, calc_fn& func){
	pair<unsigned int, const char*> retval(0,0);

	// first, let's get our data backet unpacked, please!
	stringstream ss_in;
	ss_in.write(in_buf, bufsz);
	ss_in.seekg(0);
	boost::archive::binary_iarchive ia_query(ss_in);
	mpi_query q;
	ia_query >> q;



	Result* r = func(q.calc_data->outcome, q.calc_data->data,
			q.calc_data->n_cols, q.calc_data->n_sampl, 0, q.calc_data->red_vars,
			q.class_data);

	if(r){
		r->prefix = q.calc_data->prefix;
	}

	delete q.calc_data;
	delete q.class_data;

	stringstream ss_out;
	boost::archive::binary_oarchive oa(ss_out);
	mpi_response resp;
	resp.msg_id = q.msg_id;
	resp.result = r;

	oa << resp;
	retval.first = ss_out.tellp();
	ss_out.seekg(0);
	char* buf_out = new char[retval.first];
	ss_out.read(buf_out,retval.first);
	retval.second = buf_out;

	if(r){
		delete r;
	}

	return retval;
}

}
}

