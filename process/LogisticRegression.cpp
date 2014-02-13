#include "LogisticRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <cstring>
#include <numeric>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "util/Logger.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;

using std::vector;
using std::string;
using std::set;
using std::numeric_limits;
using std::fabs;
using std::log;
using std::exp;
using std::pow;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string LogisticRegression::stepname = LogisticRegression::doRegister("logistic");

po::options_description& LogisticRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	po::options_description logreg_opts("Logistic Regression Options");

	logreg_opts.add_options()
			("odds-ratio", po::bool_switch(&show_odds), "Display odds ratios")
			("max-iterations", po::value<unsigned int>(&maxIterations)->default_value(30), "Maximum number of iterations in the logistic regression")
			;

	opts.add(logreg_opts);

	return opts;
}

void LogisticRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
}

void LogisticRegression::printVarHeader(const std::string& var_name){
	if(!show_odds){
		Regression::printVarHeader(var_name);
	}else{
		out_f << var_name << "_Pval" << sep
			  << var_name << "_OR" << sep
			  << var_name << "_SE" << sep;
	}
}

void LogisticRegression::initData(const std::string& model_str, const DataSet& ds){

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
		Utility::Logger::log_err("ERROR: Desired phenotype has only 1 unique value; logistic regression will almost certainly fail!", true);
	}

	float min_pheno = *(uniq_pheno.begin());
	float max_pheno = *(uniq_pheno.rbegin());
	float pheno_dist = max_pheno - min_pheno;

	// scale everything between 0 and 1
	for(unsigned int i=0; i<_pheno.size(); i++){
		_pheno[i] = (_pheno[i] - min_pheno) / pheno_dist;
	}

}

Regression::Result* LogisticRegression::calculate(
		const double* Y, const double* data,
		unsigned int n_cols, unsigned int n_rows, unsigned int offset,
		unsigned int n_covars){
	// Note: n_cols is the number of columns in the data vector, which is
	// 1 + # of predictor variables

	// Find the number of predictor variables in the reduced model
	Result* null_result = 0;
	unsigned int reduced_vars = n_covars + 1;

	// If this is the case, we need to find the result for running the regression
	// on the reduced model
	if(reduced_vars != 1){

		// If this is the case, we have a situation where the reduced model is
		// not quite down to our covariates, so our "new_covars" should be
		// the size of the covariates
		unsigned int new_covars = n_covars > covar_names.size() ? covar_names.size() : 0;

		// the offset is now the old offset + difference in the number of added variables
		null_result = calculate(Y, data, reduced_vars, n_rows, offset + n_cols - reduced_vars, new_covars);
	}

	Result* r = new Result();

	// val is the value of the logit function
	// deriv is the derivative of the logit
	// log_val and log_val_c are log(val) and log(1-val), respectively.
	// NOTE: the reason we do them here is to prevent underflow; for large values
	// of the exponent, we need to use approximations for log_val and log_val_c
	double val, deriv, log_val, log_val_c;

	// This is the current estimate of the parameters
	// Note: position 0 is reserved for the intercept
	r->beta_vec = new double [n_cols];
	double* beta = r->beta_vec;

	// weight vector used for IRLS procedure
	double weight[n_rows];

	gsl_matrix_const_view X = gsl_matrix_const_view_array_with_tda(data, n_rows, n_cols, offset + n_cols);

	// gsl weight vector for IRLS procedure
	gsl_vector_view w = gsl_vector_view_array(weight, n_rows);
	// gsl beta vector
	gsl_vector_view b = gsl_vector_view_array(beta, n_cols);

	// zero out the beta
	gsl_vector_set_zero (&b.vector);

	// If we have a null result, copy the beta vector for the reduced model into
	// the current beta vector - this should save a bit of time because we're
	// starting off with a better initial guess
	if(null_result){
		std::memcpy(beta, null_result->beta_vec, (n_covars + 1)*sizeof(double));
	} else {
		// Add up all the values in Y
		double sum_Y = std::accumulate(&Y[0],&Y[n_rows],0.0);
		beta[0] = log(sum_Y / (n_rows - sum_Y)); // use natural log of the ratio
	}

	// Right-hand side of the IRLS equation.  Defined to be X*w_t + S_t^-1*(y-mu_t)
	// Or, in our parlance: rhs_i = (X*beta_t)_i + 1/deriv * (y_i - val)
	gsl_vector* rhs = gsl_vector_alloc(n_rows);

	// I need these to work with the default values
	gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(n_rows, n_cols);
	gsl_matrix* cov = gsl_matrix_alloc(n_cols, n_cols);
	double tmp_chisq;



	double LLp = numeric_limits<double>::infinity(); // stores previous value of LL to check for convergence
	double LLn, LL;

	LLn = (null_result) ? null_result->log_likelihood : 0;

	unsigned int numIterations = 0;

	double TOL= 0.000000001;
	// 2*numeric_limits<double>::epsilon() ??

	while (fabs(LLp - LL) > TOL && ++numIterations < maxIterations ) {

		// First, let's initialize the RHS to X*beta_t (rhs = 1 * X * b + 0* rhs)
		gsl_blas_dgemv(CblasNoTrans, 1, &X.matrix, &b.vector, 0, rhs);

		LLp = LL;
		LL = 0;

		// add to LL for each row
		for (unsigned int i = 0; i < n_rows; i++) {

			// calculate the value of the exponent for the individual
			gsl_vector_const_view X_i = gsl_matrix_const_row(&X.matrix, i);
			double v;
			gsl_blas_ddot(&X_i.vector, &b.vector, &v);

			// At this point, v is the value of the exponent

			// max_val is the maximum value of v above that will not result in
			// loss of precision. (with a factor of 2 in there for good luck)
			static const double max_val = -log(numeric_limits<double>::epsilon());

			if(v > max_val){
				val = 1;
				deriv = exp(-v);
				// log(f(x)) obtained by Taylor series expansion of log(x) at x=1
				// note that f(x) - 1 = -exp(-x)/(1+exp(-x)) ~= -exp(-x)
				// Also, that shows the log(1-f(x)) ~= exp(-x)
				log_val = -exp(-v);
				log_val_c = -v;
			} else {
				// we won't underflow here
				val = 1/(1+exp(-v));
				log_val = log(val);
				// however, we might underflow when calculating derivatives and log
				if(-v > max_val){
					// the traditional derivative WILL underflow
					deriv = exp(v);
					log_val_c = -exp(v);
				} else{
					deriv = val*(1-val);
					log_val_c = log(1-val);
				}
			}

			// calculate LL for this ind and add to running total
			LL -= 2 *(Y[i] * log_val + (1-Y[i]) * log_val_c);

			// get the weight and update the rhs for IRLS
			weight[i] = deriv;
			gsl_vector_set(rhs, i, v + 1/deriv * (Y[i] - val));

		}

		// when this is the first iteration, set LLn (null model) to be the current value of LL
		if (null_result == 0 && numIterations == 1) {
			LLn = LL;
		}

		// Look, magic!
		// This solves the weighted least-squares problem with the weights given
		// by the derivative.  This should be a backwards-stable algorithm!
		gsl_multifit_wlinear(&X.matrix, &w.vector, rhs, &b.vector, cov, &tmp_chisq, ws);

	} // complete iteration

	r->stderr.clear();
	r->coeffs.clear();
	r->p_vals.clear();

	for(unsigned int i=1+covar_names.size(); i<n_cols; i++){
		double c = beta[i];
		double se = sqrt( gsl_matrix_get(cov, i, i));

		r->coeffs.push_back(show_odds ? exp(c) : c);
		r->stderr.push_back(se);
		// use the wald statistic to get p-values for each coefficient
		r->p_vals.push_back( gsl_cdf_chisq_Q( pow( c/se , 2) ,1) );
	}

	addResult(r, null_result);

	if (isnan(LL)) {
		r->p_val = 1.0;
		r->log_likelihood = 0.0;
	} else {
		r->p_val = gsl_cdf_chisq_Q(fabs(LLn - LL), n_cols-n_covars-1);
		r->log_likelihood = LL;
	}

	if(null_result){
		delete null_result;
	}

	gsl_vector_free(rhs);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(ws);

	return r;
}

void LogisticRegression::process(DataSet& ds){
	runRegression(ds);
}


}
}
