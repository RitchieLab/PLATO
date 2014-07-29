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
#include <cstdio>
#include <deque>
#include <utility>

#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <boost/function.hpp>

#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include "util/MPIUtils.h"

#define BOOST_IOSTREAMS_USE_DEPRECATED
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/serialization/serialization.hpp>
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "Correction.h"
#include "Encoding.h"

#include "MPIProcess.h"

#include "data/DataSet.h"

#include "util/ThreadPool.h"


namespace PLATO{

namespace Data{
class Marker;
}

namespace Analysis{

/*!
 * \brief A base class for any regression analysis - really any statistical test
 */
class Regression : virtual public MPIProcess {

public:
	class ExtraData{
	public:
		bool interactions;
		unsigned int const_covars;
		unsigned int base_covars;
		std::string sep;
		Analysis::EncodingModel encoding;
		unsigned int* extra_df_col;
		unsigned int n_extra_df_col;

		ExtraData(unsigned int n=0) : base_covars(0), extra_df_col(0), n_extra_df_col(n){
			if(n>0){
				extra_df_col = new unsigned int[n];
			}
		}

		ExtraData(const ExtraData& o){
			const_covars = o.const_covars;
			base_covars = o.base_covars;
			sep = o.sep;
			encoding = o.encoding;
			n_extra_df_col = o.n_extra_df_col;
			extra_df_col = new unsigned int[n_extra_df_col];
			std::memcpy(extra_df_col,o.extra_df_col,sizeof(unsigned int)*n_extra_df_col);
		}

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & interactions;
			ar & const_covars;
			ar & base_covars;
			ar & sep;
			ar & encoding;
			ar & n_extra_df_col;
			if(Archive::is_loading::value){
				extra_df_col = new unsigned int[n_extra_df_col];
			}

			ar & boost::serialization::make_binary_object(extra_df_col, n_extra_df_col * sizeof(unsigned int));
		}

		virtual ~ExtraData(){
			delete[] extra_df_col;
		}
	};

protected:
	/*!
	 * \brief A class that holds all data needed to process a single regression
	 */
	class Model{
	public:
		Model() : categorical(false), permute(true), id(""), n_vars(0){}
		Model(const std::vector<std::string>& tv) : traits(tv) {}

		std::string getID() const{
			if(n_vars != markers.size() + traits.size() || id.empty()){
				id = "";
				n_vars = 0;
				for(unsigned int i=0; i<markers.size(); i++){
					id += (++n_vars > 1 ? " " : "") + markers[i]->getID();
				}
				for(unsigned int i=0; i<traits.size(); i++){
					id += (++n_vars > 1 ? " " : "") + traits[i];
				}

				if(categorical){
					id += "_categ";
				}
			}

			return id;
		}

		std::vector<const PLATO::Data::Marker*> markers;
		std::vector<std::string> traits;
		bool categorical;
		bool permute;

	private:
		mutable std::string id;
		mutable unsigned int n_vars;

	};

	class ModelGenerator{
	public:

		// use this for exhaustive
		ModelGenerator(const Data::DataSet& ds) : _ds(ds) {}
		virtual ~ModelGenerator() {}

		virtual Model* next() = 0;
		virtual void reset() = 0;

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
			_mbegin(begin), _mi1(begin), _mi2(begin), _mend(end), _tbegin(trait.begin()),
			_titr(trait.begin()), _tend(trait.end()), _ti2(trait.begin()),
			_pairwise(pw), _traits(trait.size()> 0), _nomarker(nomarker) {

			if(_mi2 != _mend){
				++_mi2;
			}

			if(_ti2 != _tend){
				++_ti2;
			}
		}
		virtual ~BasicModelGenerator() {}

		virtual Model* next();
		virtual void reset();

	private:

		const M_iter _mbegin;
		M_iter _mi1;
		M_iter _mi2;
		const M_iter _mend;

		const std::set<std::string>::const_iterator _tbegin;
		std::set<std::string>::const_iterator _titr;
		const std::set<std::string>::const_iterator _tend;
		std::set<std::string>::const_iterator _ti2;

		bool _pairwise;
		bool _traits;
		bool _nomarker;
	};

	class TargetedModelGenerator : public ModelGenerator{
	public:
		// use this for targeted
		TargetedModelGenerator(const Data::DataSet& ds, const std::deque<std::string>& models) :
			ModelGenerator(ds),	_mbegin(models.begin()), _mitr(models.begin()), _mend(models.end()) {}
		virtual ~TargetedModelGenerator() {}


		virtual Model* next();
		virtual void reset();
	private:

		const std::deque<std::string>::const_iterator _mbegin;
		std::deque<std::string>::const_iterator _mitr;
		const std::deque<std::string>::const_iterator _mend;
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
		virtual void reset();

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
		Result(unsigned short n = 0) : coeffs(0), p_vals(0), stderr(0),
				submodel(0),  p_val(0),
				log_likelihood(0), r_squared(0),
				n_dropped(0), n_vars(n), converged(true) {
			if(n > 0){
				coeffs = new float[n];
				p_vals = new float[n];
				stderr = new float[n];
			}
		}
		~Result(){
			if(coeffs){delete[] coeffs;}
			if(p_vals){delete[] p_vals;}
			if(stderr){delete[] stderr;}
			if(submodel){delete submodel;}
		}

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			// get all the primitives first
			ar & p_val;
			ar & log_likelihood;
			ar & r_squared;
			ar & n_dropped;
			ar & df;
			ar & n_vars;
			ar & converged;

			// now, the STL types
			ar & prefix;
			ar & suffix;
			// NOTE: do NOT synchronize the unimodel!!
			// It will always come back empty

			// now, the arrays of data - be sure to allocate them!
			if(Archive::is_loading::value && n_vars > 0){
				coeffs = new float[n_vars];
				p_vals = new float[n_vars];
				stderr = new float[n_vars];
			}


			for(unsigned int i=0; i<n_vars; i++){
				ar & coeffs[i];
				ar & p_vals[i];
				ar & stderr[i];
			}

			ar & submodel;
			
			if(submodel && n_vars == 1){
				std::cout << "Univariate + submodel for '" << prefix << "'" << std::endl;
			}
		}

		float* coeffs;
		float* p_vals;
		float* stderr;

		Result* submodel;
		std::vector<Result*> unimodel;

		// A string to print before anything (variable IDs, MAF, etc)
		std::string prefix;
		// A string to print AFTER all of the variables, but BEFORE p-value
		std::string suffix;

		float p_val;
		float log_likelihood;
		float r_squared;

		unsigned short n_dropped;

		// the degrees of freedom in this model
		unsigned short df;

		// # of coefficients int the model
		unsigned short n_vars;

		// Did we converge (logistic regression only)?
		bool converged;

		bool operator<(const Result& o) const {return p_val < o.p_val;}

	};

	// This structure contains all the data we need to run a regression using the
	// calculate function (with the exception of the extraData)
	struct calc_matrix{
		unsigned int n_cols;
		unsigned int n_sampl;
		unsigned int red_vars;
		std::string prefix;
		double* outcome;
		double* data;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & n_cols;
			ar & n_sampl;
			ar & red_vars;
			ar & prefix;
			if(Archive::is_loading::value){
				outcome = new double[n_sampl];
				data = new double[n_sampl * n_cols];
			}
			float out_v;
			for(unsigned int i=0; i<n_sampl; i++){
				if(!Archive::is_loading::value){
					out_v = outcome[i];
				}
				ar & out_v;
				outcome[i] = out_v;
			}
			// yes, I know that data is conceptually a 2-D matrix, but
			// we're just passing pointers and data around!
			for(unsigned int i=0; i<n_sampl*n_cols; i++){
				if(!Archive::is_loading::value){
					out_v = data[i];
				}
				ar & out_v;
				data[i] = out_v;
			}
		}

		calc_matrix() : n_cols(0), n_sampl(0), red_vars(0), outcome(0), data(0)  {}

		~calc_matrix(){
			if(outcome){delete[] outcome;}
			if(data){delete[] data;}
		}

	};

public:
	/*
	 * A base class for ALL MPI messages sent
	 */
	struct mpi_data{

		virtual ~mpi_data(){}

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){}
	};


	struct mpi_permu : public mpi_data{

		unsigned long int rng_seed;
		unsigned int n_permu;
		unsigned int permu_size;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & rng_seed;
			ar & n_permu;
			ar & permu_size;
		}

	};

	struct mpi_trait : public mpi_data{

		std::vector<float> data;
		std::string description;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & data;
			ar & description;
		}

	};

	struct mpi_marker : public mpi_data{

		boost::dynamic_bitset<> data;
		std::string description;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & description;

			// serialize the data as a vector of blocks
			std::vector<boost::dynamic_bitset<>::block_type> tmp_data;
			if(!Archive::is_loading::value){
				size = data.size();
				tmp_data.reserve(data.num_blocks());
				boost::to_block_range(data,std::back_inserter(tmp_data));
			}
			ar & size;
			ar & tmp_data;

			if(Archive::is_loading::value){
				data.resize(size);
				boost::from_block_range(tmp_data.begin(), tmp_data.end(), data);
			}

		}
		
	private:
		unsigned int size;

	};

	struct mpi_weight : public mpi_data {

		//std::vector<unsigned int> _idx;
		std::vector<float> _wt;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			//ar & _idx;
			ar & _wt;
		}

	};

	struct mpi_extra : public mpi_data {

		const ExtraData* class_data;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & class_data;
		}

	};

	/*!
	 * A structure that holds everything we need to run an MPI query
	 */
	struct mpi_query : public mpi_data {
		unsigned int msg_id;
		std::vector<unsigned int> marker_idx;
		std::vector<unsigned int> trait_idx;
		bool categorical;
		bool permute;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & msg_id;
			ar & marker_idx;
			ar & trait_idx;
			ar & categorical;
			ar & permute;
		}
	};

	struct mpi_covars : public mpi_data {
		std::vector<unsigned int> covars;
		std::vector<unsigned int> const_covars;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & covars;
			ar & const_covars;
		}

	};

	struct mpi_pheno : public mpi_data {
		std::vector<float> _pheno;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
			ar & _pheno;
		}
	};

	/*
	 * A simple class indicating "prepare for broadcast"
	 */
	struct mpi_bcast : public mpi_data {
		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
		}
	};

	/*
	 * A simple struct indicating "clean up the static variables"
	 */
	struct mpi_clean : public mpi_data {
		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<mpi_data>(*this);
		}
	};

	/*
	 * An envelope class that holds a single base pointer to an mpi_data member
	 * Note that this member will be determined polymorphically based on the
	 * specific pointer being serialized
	 */
	struct mpi_envelope{
		mpi_data* msg;

		mpi_envelope() : msg(0) {}
		// please delete the data yourself!
		~mpi_envelope() {}

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & msg;
		}
	};

	/*!
	 * A structure that holds the response to a query
	 */
	struct mpi_response{
		unsigned int msg_id;
		Result* result;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & msg_id;
			ar & result;
		}

	};
public:

	Regression() :  _use_mpi(false), _lowmem(false), n_perms(0), class_data(0),  mgp(0), msg_id(0), weight_complete(true){}
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
		 * \param run_null Should we run the null model (set to false only if a permutation)
		 */
	typedef Result* (calc_fn)(const double*, const double*, unsigned int,
			unsigned int, unsigned int, unsigned int, bool, const Regression::ExtraData*);


	virtual calc_fn& getCalcFn() const = 0;

	virtual const ExtraData* getExtraData() const;

	/*virtual Result* calculate(const double* result, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars) = 0;
	*/
	float getCategoricalWeight(const PLATO::Data::Marker* m);

	virtual bool initData() = 0;
	virtual void printResults();

	virtual void printVarHeader(const std::string& var_name);
	// Use this to print extra column information, like convergence..
	virtual void printExtraHeader() {}
	virtual std::string printExtraResults(const Result& r) {return "";}

	static unsigned int findDF(const gsl_matrix* P, unsigned int reduced_vars,
			unsigned int n_dropped, unsigned int* extra_cols, unsigned int n_extra_cols);

	virtual void processResponse(unsigned int bufsz, const char* in_buf);
	virtual std::pair<unsigned int, const char*> nextQuery();

	static void calculate_MPI(unsigned int bufsz, const char* in_buf, 
		std::deque<std::pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex, 
		calc_fn& func);
		
	static void runMPIQuery(const mpi_query* mq, 
		std::deque<std::pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex,
		calc_fn& func);

private:
	Result* run(const Model& m);

	/*! \brief A function to be used for threading
	 * This function takes a ModelGenerator and iterates through it, constantly
	 * adding results to the result deque.
	 */
	void start(const Data::DataSet& ds, const std::string& outcome);

	void addResult(Result* r);
	calc_matrix* getCalcMatrix(const Model& m, const gsl_permutation* permu = 0);

	bool resetPheno(const std::string& pheno);
	void addUnivariate(Result& r, const Model& m);

	void printMarkerHeader(const std::string& var_name);
	void printHeader(unsigned int n_snp, unsigned int n_trait);

	void printResult(const Result& r, std::ostream& of);
	void printResultLine(const Result& r, std::ostream& of);

	std::pair<unsigned int, const char*> generateMsg(const Model& m);

	static void processBroadcast();
	void initMPI();
	void MPIBroadcast(std::pair<unsigned int, const char*> msg) const;
	void MPIBroadcastPheno();
	void MPIBroadcastMarker(const Data::Marker* m);
	void MPIBroadcastTrait(const std::string& t);
	void MPIBroadcastWeights() const;
	void MPIStartBroadcast() const;
	void MPIStopBroadcast() const;

public:
	typedef boost::function1<calc_matrix*, const gsl_permutation*> genPermuData;

private:

	static void initPermutations(unsigned int n_perm, unsigned int sz,
				std::deque<gsl_permutation*>& perm_list, unsigned long int permu_seed);
	static std::string getMarkerDesc(const Data::Marker* m, const std::string& sep_str);
	static void permuteData(const std::vector<float>& pheno, float*& pheno_perm,
			const std::deque<std::vector<float> >& covars, std::vector<float*>& covar_perm,
			unsigned int permuCovars, const gsl_permutation* permu);
	static unsigned int getNumCols(unsigned int n_loci, unsigned int n_trait,
			unsigned int n_covar, bool interact, bool categorical);
	static bool addDataRow(boost::multi_array_ref<double, 2>& out_data,
			const std::vector<float>& covar_vals, const std::vector<unsigned char>& geno_vals,
			const std::vector<float>& geno_weight, const std::vector<float>& trait_vals,
			const EncodingModel& enc, bool interact, bool categorical,
			unsigned int& n_samples, unsigned int& n_missing);
	static calc_matrix* getCalcMatrix(const mpi_query& mq, const gsl_permutation* permu);

	static void runPermutations(Result* r, genPermuData& perm_fn,
			const std::deque<gsl_permutation*>& permus, calc_fn& calculate, const ExtraData* ed);

private:
	//------------------------------------------------
	// General running parameters

	//! filename of output
	std::string out_fn;

	//! are we threaded?
	bool _threaded;
	//! or using MPI?
	bool _use_mpi;

	//! Do we want to reduce our memory footprint (print to file, keeping
	// only the p-values, then sort and print individual lines)
	bool _lowmem;
	boost::iostreams::stream<boost::iostreams::file_descriptor> tmp_f;

	unsigned int n_threads;
	// number of permutations to calculate p-values - set to 0 to disable
	// permutation testing.
	unsigned int n_perms;

	unsigned long int permu_seed;

	// the actual vector of permutations
	std::deque<gsl_permutation*> permutations;

	//! list of results
	std::deque<Result*> results;
	//! list of p-values (used in lowmem setting)
	std::deque<float> result_pvals;

	mutable ExtraData* class_data;

	//! Encoding scheme for the SNPs in the regression
	EncodingModel encoding;

	//------------------------------------------------
	// Data variables

	// vector of all covariates
	std::deque<std::vector<float> > _covars;

	const Data::DataSet* ds_ptr;

	//------------------------------------------------
	// Result Reporting Variables

	//! raw p-value cutoff for displaying models
	float cutoff_p;

	//! a set of correction methods to apply to the p-values
	std::set<CorrectionModel> corr_methods;

	// univariate results by marker and trait
	std::map<const PLATO::Data::Marker*, Result*> _marker_uni_result;
	std::map<std::string, Result*> _trait_uni_result;
	
	// we'll need to delete all of these univariate models
	std::deque<Result*> _uni_results;

	std::map<const PLATO::Data::Marker*, float> categ_weight;

	// mapping of column IDs to extra degrees of freedom (this comes from the weighted encoding)
	std::map<unsigned int, unsigned int> _extra_df_map;

	//------------------------------------------------
	// Model generation variables

	// Model generator to use
	ModelGenerator* mgp;
	// Iterator for the current outcome being looked at
	std::set<std::string>::const_iterator output_itr;

	//! a file of models to use
	std::vector<std::string> model_files;
	//! a list of models from the model files
	std::deque<std::string> _models;
	//! List of traits to include
	std::set<std::string> incl_traits;
	//! List of traits to exclude
	std::set<std::string> excl_traits;
	//! List of markers to include
	std::set<std::string> incl_marker_name;
	//! covariates to use in the regression
	std::set<std::string> covar_names;
	//! covariates to use in the regression DO NOT PERMUTE!
	std::set<std::string> const_covar_names;

	//! build models with traits as well?
	bool include_traits;
	//! generate one-sided pairwise models?
	bool _onesided;
	//! include interactions?
	bool interactions;
	//! exclude markers? (only use covariates?)
	bool exclude_markers;
	//! autogenerate pairwise models
	bool pairwise;
	//! Show univariate models
	bool show_uni;
	//! Do we want to do a pheWAS??
	bool _phewas;

	//-----------------------------------------------
	// Threading variables

	boost::mutex _result_mutex;
	boost::mutex _model_gen_mutex;
	boost::mutex _categ_mutex;
	boost::mutex _univar_mmutex;
	boost::mutex _univar_tmutex;

	//-----------------------------------------------
	// MPI Variables

	// The following data structures are needed for MPI to work properly
	//! A mapping of id -> models for all things currently running
	std::map<unsigned int, const Model*> work_map;
	//! mapping of results to locks that they are waiting on to add them to the
	//! result list.
	std::map<Result*, int*> post_lock_map;
	//! mapping for results->generating models for post lock results
	std::map<Result*, const Model*> post_lock_models;
	//! A mapping of currently working (or queued) model IDs to the locks they have acquired
	std::multimap<std::string, int*> work_lock_map;
	std::set<std::string> pre_lock_set;
	
	//! prefix for the outcome
	std::string mpi_prefix;

	//! a queue of models that need to be run
	std::deque<const Model*> model_queue;

	// current message ID
	unsigned int msg_id;

	// have we computed all weights?
	bool weight_complete;

	// mapping of Marker ptr to index
	std::map<const Data::Marker*, unsigned int> marker_idx_map;

	// mapping of trait to index
	std::map<const std::string, unsigned int> trait_idx_map;

	// Static variables (used only in the slaves!)
	// a deque of marker data and the "weight"
	static std::deque<std::pair<boost::dynamic_bitset<>, float> > _marker_data;
	static std::deque<std::string> _marker_desc;
	static std::deque<std::string> _trait_desc;
	static std::deque<std::vector<float> > _trait_data;
	static std::deque<std::vector<float> > _covar_data;
	static unsigned int n_const_covars;
	static std::deque<gsl_permutation*> _permu_data;
	static std::vector<float> _curr_pheno;
	static const ExtraData* _extra_data;
	static boost::shared_mutex _mpi_mutex;
	static Utility::ThreadPool _mpi_threads;

	//------------------------------------------------
	// Structures for sorting results

	struct result_sorter{
		inline bool operator() (const Result* const & x, const Result* const& y) const{
			return (x && y) ? x->p_val < y->p_val : x < y;
		}
	};

	class pval_sorter{
	public:
		pval_sorter(const std::deque<float>& v) : _v(v) {}

		bool operator() (size_t i, size_t j){
			return _v[i] < _v[j];
		}

	private:
		const std::deque<float>& _v;
	};

protected:
	//! A string of the outcome name
	std::set<std::string> outcome_names;

	//! separator to use while printing output
	std::string sep;

	// vector of all phenotypes
	std::vector<float> _pheno;

	//! output stream to print results to
	std::ofstream out_f;

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
	} else if (!_nomarker && _traits) { // Note: !_pairwise == true here
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

template <class M_iter>
void Regression::BasicModelGenerator<M_iter>::reset() {
	_mi1 = _mbegin;
	_mi2 = _mbegin;
	_titr = _tbegin;
	_ti2 = _tbegin;

	if(_mi2 != _mend){
		++_mi2;
	}

	if(_ti2 != _tend){
		++_ti2;
	}
}

}
}

#endif /* REGRESSION_H_ */
