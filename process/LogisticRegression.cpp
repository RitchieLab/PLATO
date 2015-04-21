#include "LogisticRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <cstring>
#include <numeric>

#include <boost/lexical_cast.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "util/Logger.h"
#include "util/GSLUtils.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;
using PLATO::Analysis::Encoding;
using PLATO::Utility::Logger;

using std::vector;
using std::string;
using std::set;
using std::numeric_limits;
using std::fabs;
using std::log;
using std::exp;
using std::pow;
using std::pair;
using std::map;
using std::deque;

using boost::array;

namespace po=boost::program_options;


BOOST_CLASS_EXPORT_IMPLEMENT(PLATO::ProcessLib::LogisticRegression::LogisticData)

namespace PLATO{
namespace ProcessLib{

const std::string LogisticRegression::stepname = LogisticRegression::doRegister("logistic");
const std::string LogisticRegression::MPIname = LogisticRegression::registerMPI("logistic");

po::options_description& LogisticRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	opts.add(getExtraOptions());

	return opts;
}

po::options_description LogisticRegression::getExtraOptions(){
	po::options_description logreg_opts("Logistic Regression Options");

	logreg_opts.add_options()
		("odds-ratio", po::bool_switch(&show_odds), "Display odds ratios")
		("max-iterations", po::value<unsigned int>(&maxIterations)->default_value(30), "Maximum number of iterations in the logistic regression")
		;

	return logreg_opts;
}

void LogisticRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
}

void LogisticRegression::printVarHeader(const std::string& var_name, std::ofstream& of) const{
	if(!show_odds){
		Regression::printVarHeader(var_name, of);
	}else{
		of << var_name << "_Pval" << sep
			  << var_name << "_OR" << sep
			  << var_name << "_SE" << sep;
	}
}

bool LogisticRegression::initData(){

	bool good_pheno = true;
	//transform the phenotype from [0,1]
	set<float> uniq_pheno;
	for(unsigned int i=0; i<_pheno.size(); i++){
		if(!std::isnan(_pheno[i])){
			uniq_pheno.insert(_pheno[i]);
		}
	}

	if(uniq_pheno.size() > 2){
		Utility::Logger::log_err("WARNING: Desired phenotype has more than 2 unique values; logistic regression may not be appropriate");
	}

	if(uniq_pheno.size() < 2){
		Utility::Logger::log_err("ERROR: Desired phenotype has only " +
				boost::lexical_cast<string>(uniq_pheno.size()) +
				" unique value(s); logistic regression will almost certainly fail!", outcome_names.size() <= 1);
		good_pheno = false;
	}

	float min_pheno = *(uniq_pheno.begin());
	float max_pheno = *(uniq_pheno.rbegin());
	float pheno_dist = max_pheno - min_pheno;

	// scale everything between 0 and 1
	for(unsigned int i=0; i<_pheno.size(); i++){
		_pheno[i] = (_pheno[i] - min_pheno) / pheno_dist;
	}

	return good_pheno;
}

Regression::calc_fn& LogisticRegression::getCalcFn() const{
	return (LogisticRegression::calculate);
}

const Regression::ExtraData* LogisticRegression::getExtraData() const{
	if(!class_data){
		class_data = new LogisticData(*Regression::getExtraData());
		class_data->show_odds = show_odds;
		class_data->maxIterations = maxIterations;
	}

	return class_data;
}

Regression::Result* LogisticRegression::calculate(
		const double* Y, const double* data,
		unsigned int n_cols, unsigned int n_rows, unsigned int offset,
		unsigned int n_covars, bool run_null, const Regression::ExtraData* other_data){

	const LogisticData* extra_data = dynamic_cast<const LogisticData*>(other_data);
	if(!extra_data){
		// something went VERY wrong here!
		return 0;
	}

	// Note: n_cols is the number of columns in the data vector, which is
	// 1 + # of predictor variables

	// Find the number of predictor variables in the reduced model

	unsigned int reduced_vars = n_covars + 1;
	Result* r = new Result(n_cols - 1 - extra_data->base_covars);

	// If this is the case, we need to find the result for running the regression
	// on the reduced model
	if(reduced_vars != 1){

		// If this is the case, we have a situation where the reduced model is
		// not quite down to our covariates, so our "new_covars" should be
		// the size of the covariates
		unsigned int new_covars = n_covars > extra_data->base_covars ? extra_data->base_covars : 0;

		// the offset is now the old offset + difference in the number of added variables
		if(run_null || offset == 0){
			r->submodel = calculate(Y, data, reduced_vars, n_rows,
					offset + n_cols - reduced_vars, new_covars, run_null, other_data);
		}
	}


	// val is the value of the logit function
	// deriv is the derivative of the logit
	// log_val and log_val_c are log(val) and log(1-val), respectively.
	// NOTE: the reason we do them here is to prevent underflow; for large values
	// of the exponent, we need to use approximations for log_val and log_val_c
	//double val, deriv, log_val, log_val_c;

	// This is the current estimate of the parameters
	// Note: position 0 is reserved for the intercept

	gsl_vector* beta = gsl_vector_calloc(n_cols);

	gsl_vector* weight = gsl_vector_calloc(n_rows);

	// weight vector used for IRLS procedure
	//double weight[n_rows];

	gsl_matrix_const_view data_mat = gsl_matrix_const_view_array_with_tda(data, n_rows, n_cols, offset + n_cols);

	gsl_matrix* P = gsl_matrix_alloc(n_cols, n_cols);
	// First, let's check for colinearity!
	r->n_dropped = Utility::GSLUtils::checkColinear(&data_mat.matrix, P);

	unsigned int n_indep = n_cols - r->n_dropped;

	// gsl weight vector for IRLS procedure
	//gsl_vector_view w = gsl_vector_view_array(weight, n_rows);

	// Right-hand side of the IRLS equation.  Defined to be X*w_t + S_t^-1*(y-mu_t)
	// Or, in our parlance: rhs_i = (X*beta_t)_i + 1/deriv * (y_i - val)
	gsl_vector* rhs = gsl_vector_alloc(n_rows);

	// I need these to work with the default values
	gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(n_rows, n_indep);
	gsl_matrix* cov_mat = gsl_matrix_calloc(n_cols, n_cols);
	gsl_matrix_set_all(cov_mat, -1);
	gsl_matrix_view cov_view = gsl_matrix_submatrix(cov_mat, 0, 0, n_indep, n_indep);
	gsl_matrix* cov = &cov_view.matrix;
	gsl_matrix* A = gsl_matrix_calloc(n_rows, n_cols);
	double tmp_chisq;

	// Let's perform our permutation ans set A = data * P
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &data_mat.matrix, P, 0.0, A);

	// this is the previous beta vector, for checking convergence
	gsl_vector* b_prev = gsl_vector_alloc(n_indep);
	// make this b_prev nonzero to begin
	gsl_vector_set_all(b_prev, 1.0);
	gsl_vector_view b = gsl_vector_subvector(beta, 0, n_indep);
	gsl_matrix_const_view X = gsl_matrix_const_submatrix(A, 0, 0, n_rows, n_indep);

	double LLp = numeric_limits<double>::infinity(); // stores previous value of LL to check for convergence
	double LLn, LL = 0;

	LLn = 0;

	// Add up all the values in Y
	double sum_Y = std::accumulate(&Y[0],&Y[n_rows],0.0);
	r->prefix = boost::lexical_cast<string>(sum_Y) + extra_data->sep;

	// get the exponent (and derivative, log and 1-log values) for the
	// null model (i.e., best fit of the intercept parameter)
	array<double, 4> null_v = linkFunction(log(sum_Y / (n_rows - sum_Y)));

	for(unsigned int i=0; i<n_rows; i++){
		LLn -= 2 *(Y[i] * null_v[2] + (1-Y[i]) * null_v[3]);
	}

	unsigned int numIterations = 0;

	// set the tolerances in single precision, but do work in double precision!
	double TOL= numeric_limits<float>::epsilon();
	double MAX_vec = numeric_limits<float>::max();
	// 2*numeric_limits<double>::epsilon() ??

	// check for the following:
	// 1) convergence of the likelihood
	// 1a) convergence of the beta vector (distance stored in b_prev)
	// 2) divergence of one (or more) coefficients
	// 3) maximum number of iterations
	while ((fabs(LLp - LL) > TOL*LLn || gsl_blas_dasum(b_prev) > TOL*n_cols*gsl_blas_dasum(&b.vector)) &&
		   gsl_blas_dasum(&b.vector) < MAX_vec &&
		   ++numIterations < extra_data->maxIterations ) {

		// save the old beta vector in b_prev
		gsl_vector_memcpy(b_prev, &b.vector);

		// First, let's initialize the RHS to X*beta_t (rhs = 1 * X * b + 0* rhs)
		gsl_blas_dgemv(CblasNoTrans, 1, &X.matrix, &b.vector, 0, rhs);

		LLp = LL;
		LL = 0;

		// add to LL for each row
		for (unsigned int i = 0; i < n_rows; i++) {


			// calculate the value of the exponent for the individual
			double v = gsl_vector_get(rhs, i);

			// At this point, v is the value of the exponent

			array<double, 4> v_arr = linkFunction(v);

			// calculate LL for this ind and add to running total
			LL -= 2 *(Y[i] * v_arr[2] + (1-Y[i]) * v_arr[3]);

			// get the weight and update the rhs for IRLS
			//weight[i] = v_arr[1];
			gsl_vector_set(weight, i, v_arr[1] );
			gsl_vector_set(rhs, i, v + 1/v_arr[1] * (Y[i] - v_arr[0]));

		}

		// Look, magic!
		gsl_multifit_wlinear(&X.matrix, weight, rhs, &b.vector, cov, &tmp_chisq, ws);

		// check for NaNs here
		if(std::isfinite(gsl_blas_dasum(&b.vector))){
			// get the difference between the old beta and the new beta
			gsl_vector_sub(b_prev, &b.vector);
		} else {
			// terminate the iteration, giving us the previous beta
			gsl_vector_memcpy(&b.vector, b_prev);
			// and set the "difference" to 0
			gsl_vector_set_zero(b_prev);
			// set loglikelihood to NaN
			LL = std::numeric_limits<double>::quiet_NaN();
		}

		LL += 0;

	} // complete iteration

	r->prefix += boost::lexical_cast<string>(numIterations) + extra_data->sep;


	// nonconvergence happens if:
	// -Log likelihood is not finite (inf or NaN)
	// too many iteratons
	// The current log likelihood is less than the null model
	// a submodel did not converge
	if(!std::isfinite(LL) ||
	   numIterations >= extra_data->maxIterations ||
	   LL-LLn > 0 ||
	   (r->submodel && !r->submodel->converged)){
		r->converged = false;
	}

	// OK, now time to unpermute everything!
	// Note: to unpermute, multiply by P transpose!
	// Also, we need to unpermute both the rows AND columns of cov_mat
	//b = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_matrix* _cov_work = gsl_matrix_calloc(n_cols, n_cols);
	// permute columns
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, cov_mat, P, 0.0, _cov_work);
	// permute rows
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, _cov_work, 0.0, cov_mat);
	gsl_matrix_free(_cov_work);

	gsl_vector* _bv_work = gsl_vector_alloc(n_cols);
	gsl_vector_memcpy(_bv_work, beta);
	gsl_blas_dgemv(CblasTrans, 1.0, P, _bv_work, 0.0, beta);

	// create a set of all of the removed indices
	// I'm going to re-use _bv_work from earlier to save a few bytes of memory
	gsl_vector_free(_bv_work);

	unsigned int idx_offset = 1 + extra_data->base_covars;
	for (unsigned int i = 0; i < n_cols - idx_offset; i++) {
		double c = gsl_vector_get(beta, i + idx_offset);
		double se = gsl_matrix_get(cov_mat, i + idx_offset, i + idx_offset);

		if (se >= 0) {

			r->coeffs[i] = extra_data->show_odds ? exp(c) : c;
			r->stderr[i] = sqrt(se);

			// use the wald statistic to get p-values for each coefficient
			r->p_vals[i] = gsl_cdf_chisq_Q( pow( c / r->stderr[i] , 2) ,1 + (extra_data->encoding == Encoding::WEIGHTED)) * PVAL_OFFSET;
		} else {
			// If this is true, this column was dropped from analysis!
			r->coeffs[i] = std::numeric_limits<float>::quiet_NaN();
			r->stderr[i] = std::numeric_limits<float>::quiet_NaN();
			r->p_vals[i] = std::numeric_limits<float>::quiet_NaN();
		}

	}

	r->df = findDF(P, reduced_vars, r->n_dropped, extra_data->extra_df_col, extra_data->n_extra_df_col);
	r->log_likelihood = std::isfinite(LL) ? LL : -std::numeric_limits<float>::infinity();

	Result* curr_res = r;

	string extraSuff = "";

	unsigned int df = 0;
	while (curr_res) {
		df += curr_res->df;
		pair<float, float> pv_rsq = calcPVal(r, curr_res, df, LLn);

		if (curr_res == r) {
			r->p_val = r->converged ? pv_rsq.first : (1.0f * PVAL_OFFSET);
			r->r_squared = pv_rsq.second;
		} else if (curr_res->n_vars > 0) {
			extraSuff = boost::lexical_cast<string>(pv_rsq.first * PVAL_OFFSET_RECIP) + extra_data->sep;
			break;
		}

		curr_res = curr_res->submodel;
	}


	r->suffix += extraSuff;

	gsl_vector_free(weight);
	gsl_vector_free(beta);
	gsl_vector_free(b_prev);
	gsl_vector_free(rhs);
	gsl_matrix_free(cov_mat);
	gsl_matrix_free(P);
	gsl_multifit_linear_free(ws);
	gsl_matrix_free(A);

	return r;
}

void LogisticRegression::process(DataSet& ds){
	runRegression(ds);
}

void LogisticRegression::printExtraHeader(std::ofstream& of){
	of << "Num_Cases" << sep << "N_Iter" << sep << "Converged" << sep << "Raw_LRT_pval" << sep;
}

string LogisticRegression::printExtraResults(const Result& r){
	static bool warned=false;
	if(!warned && !r.converged){
		Logger::log_err("WARNING: One or more logistic regression models did not converge");
		warned  = true;
	}
	string pv_str = boost::lexical_cast<string>(r.p_val * PVAL_OFFSET_RECIP);
	if(!r.converged){
		pair<float, float> pv_pair = calcPVal(&r, &r, r.df, 0);
		pv_str = boost::lexical_cast<string>(pv_pair.first * PVAL_OFFSET_RECIP);
	}
	return boost::lexical_cast<string>(r.converged) + sep + pv_str + sep;


}

array<double, 4> LogisticRegression::linkFunction(double v){
	// returns val, deriv (val * 1-val), log(val), 1-log(val)

	array<double, 4> retval;
	// max_val is the maximum value of v above that will not result in
	// loss of precision. (with a factor of 2 in there for good luck)
	static const double max_val = -log(numeric_limits<double>::epsilon());

	if (v > max_val) {
		retval[0] = 1;
		retval[1] = exp(-v);
		// log(f(x)) obtained by Taylor series expansion of log(x) at x=1
		// note that f(x) - 1 = -exp(-x)/(1+exp(-x)) ~= -exp(-x)
		// Also, that shows the log(1-f(x)) ~= exp(-x)
		retval[2] = -exp(-v);
		retval[3] = -v;
	} else {
		// we won't underflow here
		retval[0] = 1 / (1 + exp(-v));
		retval[2] = log(retval[0]);
		// however, we might underflow when calculating derivatives and log
		if (-v > max_val) {
			// the traditional derivative WILL underflow
			retval[1] = exp(v);
			retval[3] = -exp(v);
		} else {
			retval[1] = retval[0] * (1 - retval[0]);
			retval[3] = log(1 - retval[0]);
		}
	}
	return retval;
}

pair<float, float> LogisticRegression::calcPVal(const Result* r, const Result* curr_res, unsigned int df, float null_ll){
	pair<float, float> pv_rsq(PVAL_OFFSET,1);
	double LLn = null_ll;

	if(curr_res->submodel){
		LLn = curr_res->submodel->log_likelihood;
	}

	if (std::isfinite(r->log_likelihood) && r->log_likelihood-LLn <= 0){
		if(df != 0){
			pv_rsq.first = gsl_cdf_chisq_Q(fabs(LLn - r->log_likelihood), df) * PVAL_OFFSET;
		}
	}

	return pv_rsq;

}

void LogisticRegression::calculate_MPI(unsigned int bufsz, const char* buf,
	deque<pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex, boost::condition_variable& cv){
	Regression::calculate_MPI(bufsz, buf, result_queue, result_mutex, cv, LogisticRegression::calculate);
}

}
}
