#include "LinearRegression.h"

#include <limits>
#include <utility>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

#include "util/Sample.h"
#include "util/Marker.h"
#include "util/Logger.h"

using Methods::DataSet;
using Methods::Sample;
using Methods::Marker;
using Methods::Analysis::Regression;
using Methods::Analysis::EncodingModel;
using Methods::Analysis::CorrectionModel;

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

Regression::Result* LinearRegression::calculate(double* data, unsigned int n_cols, unsigned int n_rows){
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

	// mean and standard deviation of the coefficients (not incl. intercept!)
	vector<double> xM(n_cols -1, 0.0);
	vector<double> xSD(n_cols - 1, 0.0);

	// store xM and xSD for mean and standard deviation calculations
	for (unsigned int i = 0; i < n_rows; i++) {
		for (unsigned int j = 0; j <= n_cols-1; j++) {
			double x = data[i*n_cols + j+1];
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

	r->coeffs.clear();
	r->stderr.clear();
	r->p_vals.clear();
	// add all the coeffecients that we calculated:
	for(unsigned int i=1; i<n_cols; i++){
		r->coeffs.push_back(gsl_vector_get(beta, i));
		r->stderr.push_back(sqrt(gsl_matrix_get(cov, i, i)));
		// t-val = | beta / stderr |
		r->p_vals.push_back(2*(1-gsl_cdf_tdist_P(fabs(r->coeffs[i] / r->stderr[i]),n_rows-n_cols+1)));
	}

	double R_squared = 1 - chisq / gsl_stats_tss(Y, 1, n_rows);
	double F=(R_squared * (n_rows - n_cols))/((1-R_squared) * (n_cols - 1));

	r->p_val = 1-gsl_cdf_fdist_P(F,n_cols-1,n_rows-n_cols);

	// calulate the sum of the squares of the residuals
	double sum_sq_resid = gsl_stats_tss_m(resid->data, 1, n_rows, 0);
	r->log_likelihood = 0.5 * (-n_rows * (log(2*M_PI)+1 - log(n_rows) + log(sum_sq_resid)));

	gsl_vector_free(beta);
	gsl_vector_free(resid);
	gsl_matrix_free(cov);

	return r;
}

void LinearRegression::process(DataSet& ds){
	runRegression(ds);
}


}
