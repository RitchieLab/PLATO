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

Regression::Result* LinearRegression::calculate(double* data, unsigned int n_cols, unsigned int n_rows, const Result* null_result){
	// Note: n_cols is the number of columns in the data vector, which is
	// 1 + # of predictor variables

	Result* r = new Result();

	// We want to calculate the best fit for X*b = y
	double Y[n_rows];

	gsl_matrix_view X = gsl_matrix_view_array (data, n_rows, n_cols);
	gsl_vector_view y_vec = gsl_vector_view_array(Y, n_rows);
	gsl_vector* beta = gsl_vector_alloc(n_cols);
	gsl_vector* resid = gsl_vector_alloc(n_rows);
	gsl_matrix* cov = gsl_matrix_alloc(n_cols, n_cols);
	double chisq;

	for(unsigned int i=0; i<n_rows; i++){
		Y[i] = data[i*n_cols];
		data[i*n_cols] = 1;
	}

	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(n_rows, n_cols);

	gsl_multifit_linear(&X.matrix, &y_vec.vector, beta, cov, &chisq, ws);
	gsl_multifit_linear_residuals(&X.matrix, &y_vec.vector, beta, resid);

	// OK, now we have our coefficients!  It's that easy!

	r->coeffs.clear();
	r->stderr.clear();
	r->p_vals.clear();
	// add all the coeffecients that we calculated:
	for(unsigned int i=1; i<n_cols; i++){
		r->coeffs.push_back(gsl_vector_get(beta, i));
		r->stderr.push_back(sqrt(gsl_matrix_get(cov, i, i)));
		// t-val = | beta / stderr |
		r->p_vals.push_back(2*(1-gsl_cdf_tdist_P(fabs(r->coeffs[i-1] / r->stderr[i-1]),n_rows-n_cols+1)));
	}

	double tss = gsl_stats_tss(Y, 1, n_rows);
	double null_rss = tss;
	int n_null = 0;
	if(null_result != 0){
		null_rss = tss*(1-null_result->r_squared);
		n_null = null_result->coeffs.size();
	}


	r->r_squared = 1 - chisq / tss;
	double F=((null_rss - chisq) * (n_rows - n_cols))/(chisq * (n_cols - 1 - n_null));

	r->p_val = 1-gsl_cdf_fdist_P(F,n_cols-1-n_null,n_rows-n_cols);

	// calulate the sum of the squares of the residuals
	//double sum_sq_resid = gsl_stats_tss_m(resid->data, 1, n_rows, 0);
	double ll = 0.5 * (-n_rows * (log(2*M_PI)+1 - log(n_rows) + log(chisq)));
	if(null_result){
		r->log_likelihood = -2 * (null_result->log_likelihood - ll);
	}else{
		r->log_likelihood = ll;
	}

	double pv_test = gsl_cdf_chisq_Q(r->log_likelihood,1);

	//r->log_likelihood = r_full->log_likelihood;
	//r->p_val = r_full->p_val;
	//r->p_val = gsl_cdf_chisq_Q(r->log_likelihood,1);

	gsl_vector_free(beta);
	gsl_vector_free(resid);
	gsl_matrix_free(cov);

	return r;
}

void LinearRegression::process(DataSet& ds){
	runRegression(ds);
}


}
}
