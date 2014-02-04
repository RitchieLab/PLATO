#include "LogisticRegression.h"

#include <limits>
#include <utility>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "util/Sample.h"
#include "util/Marker.h"
#include "util/Logger.h"

using Methods::DataSet;
using Methods::Sample;
using Methods::Marker;
using Methods::Analysis::Regression;
using Methods::Analysis::EncodingModel;

using std::vector;
using std::deque;
using std::string;
using std::set;
using std::numeric_limits;
using std::pair;
using std::fabs;
using std::log;
using std::exp;

namespace po=boost::program_options;

namespace ProcessLib{

const std::string LogisticRegression::stepname = LogisticRegression::doRegister("logistic");

po::options_description& LogisticRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	po::options_description logreg_opts("Logistic Regression Options");

	logreg_opts.add_options()
			//("odds-ratio", po::bool_switch(&show_odds), "Display odds ratios")
			("max-iterations", po::value<unsigned int>(&maxIterations)->default_value(30), "Maximum number of iterations in the logistic regression")
			;

	opts.add(logreg_opts);

	return opts;
}

void LogisticRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
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

Regression::Result* LogisticRegression::calculate(double* data, unsigned int n_cols, unsigned int n_rows){
	// Note: n_cols is the number of columns in the data vector, which is
	// 1 + # of predictor variables

	Result* r = new Result();


	double x;

	// val is the value of the logit function
	// deriv is the derivative of the logit
	// log_val and log_val_c are log(val) and log(1-val), respectively.
	// NOTE: the reason we do them here is to prevent underflow; for large values
	// of the exponent, we need to use approximations for log_val and log_val_c
	double val, deriv, log_val, log_val_c;

	// mean and standard deviation of the coefficients (not incl. intercept!)
	vector<double> xM(n_cols -1, 0.0);
	vector<double> xSD(n_cols - 1, 0.0);

	// This is the current estimate of the parameters
	// Note: position 0 is reserved for the intercept
	double beta[n_cols];

	// weight vector used for IRLS procedure
	double weight[n_rows];
	double Y[n_rows];

	gsl_matrix_view X = gsl_matrix_view_array(data, n_cols, n_rows);

	// gsl weight vector for IRLS procedure
	gsl_vector_view w = gsl_vector_view_array(weight, n_rows);
	// gsl beta vector
	gsl_vector_view b = gsl_vector_view_array(beta, n_cols);

	// Right-hand side of the IRLS equation.  Defined to be X*w_t + S_t^-1*(y-mu_t)
	// Or, in our parlance: rhs_i = (X*beta_t)_i + 1/deriv * (y_i - val)
	gsl_vector* rhs = gsl_vector_alloc(n_rows);

	// I need these to work with the default values
	gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(n_rows, n_cols);
	gsl_matrix* cov = gsl_matrix_alloc(n_cols, n_cols);
	double tmp_chisq;

	// store xM and xSD for mean and standard deviation calculations
	for (unsigned int i = 0; i < n_rows; i++) {

		for (unsigned int j = 0; j <= n_cols-1; j++) {
			x = data[i*n_cols + j+1];
			xM[j] += x;
			xSD[j] += x*x;
		}

		Y[i] = data[i*n_cols];
		data[i*n_cols] = 1;
	}

	// calculate mean and standard deviation
	for (unsigned int j = 0; j <= n_cols-1; j++) {
		xM[j] /= n_rows;
		xSD[j] /= n_rows;
		xSD[j] = sqrt(fabs(xSD[j] - xM[j] * xM[j]));
	}

	// adjusts X values using the mean and standard deviation values
	/*
	for (unsigned int i = 0; i < n_rows; i++) {
		for (unsigned int j = 0; j <= n_cols-1; j++) {
			// This subtracts the mean, then divides by the stddev
			(data[i*n_cols + j+1] -= xM[j]) /= xSD[j] ;
		}
	}
	*/

	double sY1 = 0;
	for(unsigned int i=0; i< n_rows; i++){
		sY1 += Y[i];
	}
	double sY0 = n_rows - sY1;

	beta[0] = log(sY1 / sY0); // use natural log of the ratio

	double LLp = numeric_limits<double>::infinity(); // stores previous value of LL to check for convergence
	double LLn = 0, LL = 0;
	unsigned int numIterations = 0;

	while (fabs(LLp - LL) > 0.0000001 && ++numIterations > maxIterations ) {

		// First, let's initialize the RHS to X*beta_t (rhs = 1 * X * b + 0* rhs)
		gsl_blas_dgemv(CblasNoTrans, 1, &X.matrix, &b.vector, 0, rhs);

		LLp = LL;
		LL = 0;

		// add to LL for each row
		for (unsigned int i = 0; i < n_rows; i++) {

			// calculate the value of the exponent for the individual
			double v = beta[0];
			for (unsigned int j = 1; j <= n_cols; j++) {
				v += v + beta[j] * data[i*n_cols + j];
			}

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
			gsl_vector_set(rhs, i, gsl_vector_get(rhs, i) + 1/deriv * (Y[i] - val));

		}

		// when this is the first iteration, set LLn (null model) to be the current value of LL
		if (numIterations == 1) {
			LLn = LL;
		}

		// Look, magic!
		// This solves the weighted least-squares problem with the weights given
		// by the derivative.  This should be a backwards-stable algorithm!
		gsl_multifit_wlinear(&X.matrix, &w.vector, rhs, &b.vector, cov, &tmp_chisq, ws);

	} // complete iteration

	// calculate p values for the coefficients
	// interaction coefficient for all loci is the first one
	/*
	for (unsigned int j = 1; j <= n_cols; j++) {
		beta[j] /= xSD[j-1];
		beta[0] -= beta[j] * xM[j-1];
	}
	*/

	r->stderr.clear();
	r->coeffs.clear();
	r->p_vals.clear();
	for (unsigned int j = 1; j <= n_cols; j++) {
		r->coeffs.push_back(beta[j]);
		r->stderr.push_back(sqrt(gsl_matrix_get(cov, j, j)/weight[j]));
		// use the wald statistic to get p-values for each coefficient
		r->p_vals.push_back(gsl_cdf_chisq_Q(fabs(r->coeffs[j] / r->stderr[j]),1));
	}

	if (isnan(LL)) {
		r->p_val = 1.0;
		r->log_likelihood = 0.0;
	} else {
		r->p_val = gsl_cdf_chisq_Q(fabs(LLn - LL), n_cols);
		r->log_likelihood = LL - LLn;
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
