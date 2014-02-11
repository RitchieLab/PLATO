/*
 * Regression.h
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#ifndef METHODS_ANALYSIS_REGRESSION_H
#define METHODS_ANALYSIS_REGRESSION_H

#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <boost/program_options.hpp>

#include "Correction.h"
#include "Encoding.h"

#include "data/DataSet.h"

namespace PLATO{

namespace Data{
class Marker;
}

namespace Analysis{

/*!
 * \brief A base class for any regression analysis - really any statistical test
 */
class Regression {

protected:
	/*!
	 * \brief A class that holds all data needed to process a single regression
	 */
	class Model{
	public:
		Model() : interact(false), categorical(false), reduce(true){}
		Model(const std::vector<std::string>& tv) : traits(tv) {}

		std::vector<const PLATO::Data::Marker*> markers;
		std::vector<std::string> traits;

		bool interact;
		bool categorical;
		bool reduce;

	};

	class ModelGenerator{
	public:

		// use this for exhaustive
		ModelGenerator(const PLATO::Data::DataSet& ds,
				const std::set<std::string>& trait, bool pw,
				bool nomarker) :
			_ds(ds), _mi1(ds.beginMarker()), _mi2(ds.beginMarker()),
			_titr(trait.begin()), _tend(trait.end()), _ti2(trait.begin()),
			_targeted(false), _pairwise(pw), _traits(trait.size()> 0), _nomarker(nomarker) {
		}

		// use this for targeted
		ModelGenerator(const PLATO::Data::DataSet& ds, const std::deque<std::string>& models, bool interact) :
			_ds(ds), _mi1(ds.beginMarker()), _mi2(ds.beginMarker()),
			_mitr(models.begin()), _mend(models.end()),
			_targeted(true) {

		}

		// Return a new Model* object, or NULL if finished
		Model* operator() ();
	private:

		const PLATO::Data::DataSet& _ds;

		PLATO::Data::DataSet::const_marker_iterator _mi1;
		PLATO::Data::DataSet::const_marker_iterator _mi2;

		std::set<std::string>::const_iterator _titr;
		std::set<std::string>::const_iterator _tend;
		std::set<std::string>::const_iterator _ti2;

		std::deque<std::string>::const_iterator _mitr;
		std::deque<std::string>::const_iterator _mend;

		bool _targeted;
		bool _pairwise;
		bool _traits;
		bool _nomarker;

	};

	/*!
	 * \brief A class that holds the result of the regression and can be used
	 * to print the results.
	 */
	class Result{
	public:
		std::deque<float> coeffs;
		std::deque<float> p_vals;
		std::deque<float> stderr;

		float p_val;
		float log_likelihood;
		float r_squared;

		// A string to print before anything (variable IDs, MAF, etc)
		std::string prefix;
		// A string to print AFTER all of the variables, but BEFORE p-value
		std::string suffix;

		bool operator<(const Result& o) const {return p_val < o.p_val;}
	};

public:

	Regression() : out_f(NULL) {}
	virtual ~Regression();

	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);
	void parseOptions(const boost::program_options::variables_map& vm);

	//! iterate through and run the regressions we want.
	void runRegression(const PLATO::Data::DataSet& ds);

	static Model* parseModelStr(const std::string& model_str, const PLATO::Data::DataSet& ds);

protected:

	virtual Result* calculate(double* data, unsigned int n_cols, unsigned int n_rows, const Result* null_result=0) = 0;
	float getCategoricalWeight(const PLATO::Data::Marker* m, const PLATO::Data::DataSet& ds);

	virtual void initData(const std::string& model_str, const PLATO::Data::DataSet& ds) {}
	virtual void printResults();

	virtual void printVarHeader(const std::string& var_name);

private:
	virtual Result* run(const Model* m, const PLATO::Data::DataSet& ds, Result* null_model=0);


private:
	//! a file of models to use
	std::vector<std::string> model_files;
	//! a list of models from the model files
	std::deque<std::string> _models;
	//! filename of output
	std::string out_fn;

	//! List of traits to exclude
	std::set<std::string> excl_traits;

	//! build models with traits as well?
	bool include_traits;

	std::map<const PLATO::Data::Marker*, float> categ_weight;

	struct result_sorter{
		inline bool operator() (const Result* const & x, const Result* const& y) const{
			return (x && y) ? x->p_val < y->p_val : x < y;
		}
	};

protected:

	//! List of traits to include
	std::set<std::string> incl_traits;

	//! a set of correction methods to apply to the p-values
	std::set<CorrectionModel> corr_methods;
	//! covariates to use in the regression
	std::set<std::string> covar_names;
	//! A string of the outcome name
	std::string outcome_name;

	//! separator to use while printing output
	std::string sep;

	//! raw p-value cutoff for displaying models
	float cutoff_p;
	//! include interactions?
	bool interactions;
	//! exclude markers? (only use covariates?)
	bool exclude_markers;
	//! autogenerate pairwise models
	bool pairwise;
	//! Show univariate models
	bool show_uni;

	//! Encoding scheme for the SNPs in the regression
	EncodingModel encoding;

	// vector of all phenotypes
	std::vector<float> _pheno;
	// vector of all covariates
	std::vector<std::vector<float> > _covars;

	//! list of results
	std::deque<Result*> results;

	//! output stream to print results to
	std::ofstream out_f;

	// univariate results by marker and trait
	std::map<const PLATO::Data::Marker*, Result*> _marker_uni_result;
	std::map<std::string, Result*> _trait_uni_result;

};

}
}



#endif /* REGRESSION_H_ */
