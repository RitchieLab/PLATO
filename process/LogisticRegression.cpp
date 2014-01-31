#include "LogisticRegression.h"

#include <limits>
#include <utility>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

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
			("odds-ratio", po::bool_switch(&show_odds), "Display odds ratios")
			("max-iterations", po::value<unsigned int>(&maxIterations)->default_value(30), "Maximum number of iterations in the logistic regression")
			;

	opts.add(logreg_opts);

	return opts;
}

void LogisticRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
}

void LogisticRegression::initData(const std::string& model_str, const DataSet& ds){
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

Regression::Result* LogisticRegression::calculate(float* data, unsigned int n_cols, unsigned int n_rows){
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
	vector<double> beta(n_cols, 0);

	// standard errors of the coefficients
	vector<double> SEP(n_cols);

	// Use this matrix for updating parameter estimates
	// Note that this is related to the hessian
	double Hess[n_cols][n_cols];

	// We're going to use this matrix "view" to do the solve that we need
	gsl_matrix_view T = gsl_matrix_view_array (Hess[0], n_cols, n_cols);

	// This will hold the RHS of the solver
	gsl_vector* c = gsl_vector_alloc(n_cols);
	gsl_vector* beta_add = gsl_vector_alloc(n_cols);

	// Just need this for QR decomposition
	gsl_vector* tau = gsl_vector_alloc(n_cols);

	// store xM and xSD for mean and standard deviation calculations
	for (unsigned int i = 0; i < n_rows; i++) {

		for (unsigned int j = 0; j <= n_cols-1; j++) {
			x = data[i*n_cols + j+1];
			xM[j] += x;
			xSD[j] += x*x;
		}
	}

	// calculate mean and standard deviation
	for (unsigned int j = 0; j <= n_cols-1; j++) {
		xM[j] /= n_rows;
		xSD[j] /= n_rows;
		xSD[j] = sqrt(fabs(xSD[j] - xM[j] * xM[j]));
	}

	// adjusts X values using the mean and standard deviation values
	for (unsigned int i = 0; i < n_rows; i++) {
		for (unsigned int j = 0; j <= n_cols-1; j++) {
			// This subtracts the mean, then divides by the stddev
			(data[i*n_cols + j+1] -= xM[j]) /= xSD[j] ;
		}
	}

	float sY1 = 0;
	for(unsigned int i=0; i< n_rows; i++){
		sY1 += data[i*n_cols];
	}
	float sY0 = n_rows - sY1;

	beta[0] = log(sY1 / sY0); // use natural log of the ratio

	double LLp = numeric_limits<double>::infinity(); // stores previous value of LL to check for convergence
	double LLn = 0, LL = 0;
	unsigned int numIterations = 0;

	while (fabs(LLp - LL) > 0.0000001 && ++numIterations > maxIterations ) {

		LLp = LL;
		LL = 0;

		// zero out Hess for this iteration
		std::memset(Hess[0],0,n_cols * n_cols * sizeof(double));

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
			LL -= 2 *(data[i*n_cols] * log_val + (1-data[i*n_cols]) * log_val_c);

			// We're calculating X^T * W * X for the iterations here

			// First, we have the intercept terms for which X == 1
			for (unsigned int j = 0; j < n_cols; j++) {
				Hess[0][j] = (Hess[j][0] += deriv);
			}

			// Now, we have to go through our data!
			for (unsigned int j = 1; j < n_cols; j++){
				for (unsigned int k = j; k < n_cols; k++) {
					Hess[j][k] = (Hess[k][j] += data[i*n_cols + j] * deriv * data[i*n_cols + k]);
				}
			}
		}

		// when this is the first iteration, set LLn (null model) to be the current value of LL
		if (numIterations == 1) {
			LLn = LL;
		}

		// OK, I want to find b=T^-1*X'*(Y-val)
		// OR... I can solve the equation (for v): Tv = X'*(Y-val)

		// First, clear the RHS:
		gsl_vector_set_zero(c);
		// set the c = X'*(Y-val)

		// Now, do a QR factorization on the "Hess" or T matrix:
		gsl_linalg_QR_decomp(&T.matrix, tau);
		gsl_linalg_QR_solve(&T.matrix, tau, c, beta_add);

		for(unsigned int j=0; j<n_cols; j++){
			beta[j] += gsl_vector_get(beta_add, j);
		}

	} // complete iteration

	// calculate p values for the coefficients
	// interaction coefficient for all loci is the first one
	for (unsigned int j = 1; j <= n_cols; j++) {
		beta[j] /= xSD[j];
		beta[0] -= beta[j] * xM[j];
	}

	for (unsigned int j = 0; j <= n_cols; j++) {
		SEP[j] = sqrt(Hess[j][j] / xSD[j]);
	}

	r->coeffs.clear();
	// add all the coeffecients that we calculated:
	r->coeffs.insert(r->coeffs.begin(), beta.begin(), beta.end());

	r->stderr.clear();
	r->stderr.insert(r->stderr.begin(), SEP.begin(), SEP.end());

	// use the wald statistic to get p-values for each coefficient
	r->p_vals.clear();
	for (unsigned int i = 0; i < n_cols; i++) {
		r->p_vals.push_back(gsl_cdf_chisq_Q(fabs(beta[i] / SEP[i]), 1));
	}

	if (isnan(LL)) {
		r->p_val = 1.0;
		r->log_likelihood = 0.0;
	} else {
		r->p_val = gsl_cdf_chisq_Q(fabs(LLn - LL), n_cols);
		r->log_likelihood = LL - LLn;
	}


	gsl_vector_free(c);
	gsl_vector_free(beta_add);
	gsl_vector_free(tau);

	return r;
}

void LogisticRegression::printResults(){
	// Perhaps sort the deque of results based on overall p-values

	// Now, let's do some multiple test correction!


	for(unsigned int i=0; i<results.size(); i++){
		// If we've gone over the threshold, stop printing!
		if(results[i]->p_val > cutoff_p){
			break;
		}
		// print a single line
	}
}

void LogisticRegression::process(DataSet& ds){
	runRegression(ds);
}


}
