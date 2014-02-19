#include "LinearRegression.h"

#include <limits>
#include <utility>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

#include "util/Logger.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;
using PLATO::Analysis::Encoding;

using std::vector;
using std::string;
using std::fabs;
using std::log;
using std::exp;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string LinearRegression::stepname = LinearRegression::doRegister("linear");

po::options_description& LinearRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	po::options_description logreg_opts("Linear Regression Options");

	logreg_opts.add_options()
			;

	opts.add(logreg_opts);

	return opts;
}

void LinearRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
}

Regression::Result* LinearRegression::calculate(
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
		null_result = calculate(Y, data, reduced_vars, n_rows, offset + n_cols - (reduced_vars), new_covars);
	}

	Result* r = new Result();

	r->beta_vec = new double[n_cols];
	// We want to calculate the best fit for X*b = y
	//double Y[n_rows];

	gsl_matrix_const_view X = gsl_matrix_const_view_array_with_tda(data, n_rows, n_cols, offset + n_cols);
	gsl_vector_const_view y_vec = gsl_vector_const_view_array(Y, n_rows);
	gsl_vector_view bv = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_vector* beta = &bv.vector;
	gsl_vector* resid = gsl_vector_alloc(n_rows);
	gsl_matrix* cov = gsl_matrix_alloc(n_cols, n_cols);
	double chisq;

	/*for(unsigned int i=0; i<n_rows; i++){
		Y[i] = data[i*n_cols];
		data[i*n_cols] = 1;
	}*/

	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(n_rows, n_cols);

	gsl_multifit_linear(&X.matrix, &y_vec.vector, beta, cov, &chisq, ws);
	gsl_multifit_linear_residuals(&X.matrix, &y_vec.vector, beta, resid);

	// OK, now we have our coefficients!  It's that easy!

	r->coeffs.clear();
	r->stderr.clear();
	r->p_vals.clear();
	// add all the non-covariate coefficients
	for(unsigned int i=1+covar_names.size(); i<n_cols; i++){
		double c = gsl_vector_get(beta, i);
		double se = sqrt(gsl_matrix_get(cov, i, i));
		r->coeffs.push_back(c);
		r->stderr.push_back(se);
		// t-val = | beta / stderr |
		r->p_vals.push_back(2*(1-gsl_cdf_tdist_P(fabs(c / se),n_rows-n_cols+1)));
	}

	addResult(r, null_result);


	double tss = gsl_stats_tss(Y, 1, n_rows);
	r->r_squared = 1 - chisq/tss;
	double null_rss = tss;
	if (null_result){
		null_rss *= 1-null_result->r_squared;
	}

	double F=((null_rss - chisq) * (n_rows - n_cols))/(chisq * (n_cols - reduced_vars));

	// We want to see if there are extra degrees of freedom, which can happen in
	// the case of the "categorical" model
	// we have an extra df per marker in the categorically encoded model
	// We only have markers if we are not excluding markers and the
	// number of columns is at least as many as the number of covariates
	// (i.e. this isn;t the "null" model)
	unsigned int extra_df = (encoding == Encoding::CATEGORICAL)
			* (!interactions || offset != 0)
			* (!exclude_markers) * (n_cols > covar_names.size() + 1)
			* (1 + pairwise);


	r->p_val = 1-gsl_cdf_fdist_P(std::max(0.0, F),n_cols-reduced_vars+extra_df,n_rows-n_cols-extra_df);

	// I have no idea if the log_likelihood is correct!!
	r->log_likelihood = 0.5 * (-n_rows * (log(2*M_PI)+1 - log(n_rows) + log(chisq)));

	double pv_test = gsl_cdf_chisq_Q(r->log_likelihood,1);

	//r->log_likelihood = r_full->log_likelihood;
	//r->p_val = r_full->p_val;
	//r->p_val = gsl_cdf_chisq_Q(r->log_likelihood,1);

	// Make sure to clean up after yourself!
	if(null_result){
		delete null_result;
	}

	gsl_vector_free(resid);
	gsl_matrix_free(cov);

	return r;
}

void LinearRegression::process(DataSet& ds){
	runRegression(ds);
}


}
}
