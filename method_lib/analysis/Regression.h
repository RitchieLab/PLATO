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
#include <boost/thread.hpp>

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
		Model() : categorical(false){}
		Model(const std::vector<std::string>& tv) : traits(tv) {}

		std::vector<const PLATO::Data::Marker*> markers;
		std::vector<std::string> traits;

		bool categorical;

	};

	class ModelGenerator{
	public:

		// use this for exhaustive
		ModelGenerator(const Data::DataSet& ds) : _ds(ds) {}
		virtual ~ModelGenerator() {}

		virtual Model* next() = 0;

		// Return a new Model* object, or NULL if finished
		Model* operator() () {return next();};
	protected:

		const Data::DataSet& _ds;

	};

	template <class M_iter>
	class BasicModelGenerator : public ModelGenerator{

	public:
		BasicModelGenerator(const PLATO::Data::DataSet& ds,
				const M_iter& begin, const M_iter& end,
				const std::set<std::string>& trait, bool pw,
				bool nomarker) : ModelGenerator(ds),
			_mbegin(begin), _mi1(begin), _mi2(begin), _mend(end),
			_titr(trait.begin()), _tend(trait.end()), _ti2(trait.begin()),
			_pairwise(pw), _traits(trait.size()> 0), _nomarker(nomarker) {

			if(_mi2 != _mend){
				++_mi2;
			}
		}
		virtual ~BasicModelGenerator() {}

		virtual Model* next();

	private:

		M_iter _mbegin;
		M_iter _mi1;
		M_iter _mi2;
		M_iter _mend;

		std::set<std::string>::const_iterator _titr;
		std::set<std::string>::const_iterator _tend;
		std::set<std::string>::const_iterator _ti2;

		bool _pairwise;
		bool _traits;
		bool _nomarker;
	};

	class TargetedModelGenerator : public ModelGenerator{
	public:
		// use this for targeted
		TargetedModelGenerator(const Data::DataSet& ds, const std::deque<std::string>& models) :
			ModelGenerator(ds),	_mitr(models.begin()), _mend(models.end()) {}
		virtual ~TargetedModelGenerator() {}


		virtual Model* next();
	private:

		std::deque<std::string>::const_iterator _mitr;
		std::deque<std::string>::const_iterator _mend;
	};

	class OneSidedModelGenerator : public ModelGenerator{

	public:
		OneSidedModelGenerator(const Data::DataSet& ds,
				const std::set<const Data::Marker*>& m_set,
				const std::set<std::string>& trait,
				const std::set<std::string>& all_traits,
				bool nomarker ) : ModelGenerator(ds), _m_set(m_set),
			_t_set(trait), _tall_set(all_traits),
			_mi1(m_set.begin()), _mi2(ds.beginMarker()),
			_titr(trait.begin()), _ti2(all_traits.begin()),
			_traits(trait.size() > 0),	_nomarker(nomarker){}
		virtual ~OneSidedModelGenerator() {}

		virtual Model* next();

		private:
			const std::set<const Data::Marker*>& _m_set;
			const std::set<std::string>& _t_set;
			const std::set<std::string>& _tall_set;

			std::set<const Data::Marker*>::const_iterator _mi1;
			Data::DataSet::const_marker_iterator _mi2;

			std::set<std::string>::const_iterator _titr;
			std::set<std::string>::const_iterator _ti2;

			std::set<std::string> _t_processed;
			std::set<const Data::Marker*> _m_processed;

			// NOTE: we ALWAYS want pairwise with this generator!
			bool _traits;
			bool _nomarker;
	};

	/*!
	 * \brief A class that holds the result of the regression and can be used
	 * to print the results.
	 */
	class Result{
	public:
		Result() : beta_vec(0), n_dropped(0), converged(true) {}
		~Result(){ if(beta_vec){delete[] beta_vec;}}

		std::deque<float> coeffs;
		std::deque<float> p_vals;
		std::deque<float> stderr;

		float p_val;
		float log_likelihood;
		float r_squared;

		// An array of beta values (including intercept + covariate)
		// we want to use this as a starting point for iterations of expanded
		// models
		double* beta_vec;

		// A string to print before anything (variable IDs, MAF, etc)
		std::string prefix;
		// A string to print AFTER all of the variables, but BEFORE p-value
		std::string suffix;

		unsigned int n_dropped;

		// Did we converge (logistic regression only)?
		bool converged;

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

	/*! \brief Runs the regression calculation
	 * This function actually runs the regression (or appropriate statistical
	 * test) calculation based on the raw data passed in.
	 * \param result A n_rows length vector containing the output vector
	 * \param data A n_rows x n_cols matrix containing the independent values
	 * for the regression.  Note that the data should be in the following order:
	 * - intercept (a column of all 1's)
	 * - covariates
	 * - main effects
	 * - interactions
	 * Also note that the data is arranged linearly in memory, with the following
	 * relationship: data[i][j] = data[i*(n_cols + offset) + j]
	 * \param n_rows The number of rows in the dataset
	 * \param n_cols The number of independent variables in the dataset
	 * \param offset The number of extra columns in the data not included in the
	 * particular model that we are testing.
	 * \param n_covars The number of covariates in the "reduced" model to test
	 * significance against.  Note that if n_covars > covar_names.size() + 1, then
	 * this method will call calculate again.  For an example, see the definition
	 * in LinearRegression.
	 */
	virtual Result* calculate(const double* result, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars=0) = 0;
	float getCategoricalWeight(const PLATO::Data::Marker* m, const PLATO::Data::DataSet& ds);

	virtual void initData(const std::string& model_str, const PLATO::Data::DataSet& ds) {}
	virtual void printResults();

	virtual void printVarHeader(const std::string& var_name);
	// Use this to print extra column information, like convergence..
	virtual void printExtraHeader() {}
	virtual void printExtraResults(const Result& r) {}

	void addResult(Result* curr_result, const Result* null_result);

private:
	Result* run(const Model* m, const Data::DataSet& ds);

	/*! \brief A function to be used for threading
	 * This function takes a ModelGenerator and iterates through it, constantly
	 * adding results to the result deque.
	 */
	void start(ModelGenerator& mg, const Data::DataSet& ds);

	void printMarkerHeader(const std::string& var_name);
	void printHeader(unsigned int n_snp, unsigned int n_trait);

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

	//! generate one-sided pairwise models?
	bool _onesided;

	//! are we threaded?
	bool _threaded;

	//! List of markers to include
	std::set<std::string> incl_marker_name;

	unsigned int n_threads;

	boost::mutex _result_mutex;
	boost::mutex _model_gen_mutex;

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

template <class M_iter>
Regression::Model* Regression::BasicModelGenerator<M_iter>::next() {

	Model* m = 0;

	// We must want exhaustive pairwise models
	if (_pairwise) {

		// We want Env vars
		if (_traits) {
			// we want exhaustive EnvxEnv
			if (_nomarker) {
				if (++_ti2 == _tend && _titr != _tend && ++_titr != _tend) {
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

				if (_mi2 == _mend) {
					if (_mi1 != _mend && ++_mi1 != _mend) {
						_mi2 = _mi1;
						if (++_mi2 ==_mend) {
							_mi1 = _mbegin;
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
		} else {

			if (_mi2 == _mend) {
				if (_mi1 != _mend && ++_mi1 != _mend){
					_mi2 = _mi1;
					++_mi2;
				}
			}
			if (_mi2 != _mend) {
				m = new Model();
				m->markers.push_back(*_mi1);
				m->markers.push_back(*_mi2);
				++_mi2;
			}

		}
		// We want Marker x Trait models (i.e. GxE)
	} else if (_traits) { // Note: !_pairwise == true here
		if (_mi1 == _mend) {
			_mi1 = _mbegin;
			++_titr;
		}
		if (_mi1 != _mend && _titr != _tend) {
			m = new Model();
			m->markers.push_back(*_mi1);
			m->traits.push_back(*_titr);
			++_mi1;
		}
		// Must want single variable models
	} else {
		// Env only
		if (_nomarker) {
			if (_titr != _tend) {
				m = new Model();
				m->traits.push_back(*_titr);
				++_titr;
			}
			// SNP only
		} else {
			if (_mi1 != _mend) {
				m = new Model();
				m->markers.push_back(*_mi1);
				++_mi1;
			}
		}
	}

	// NOTE: m == 0 if no more models can be generated!
	return m;

}

}
}



#endif /* REGRESSION_H_ */
