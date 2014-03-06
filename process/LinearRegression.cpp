#include "LinearRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <set>
#include <algorithm>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "util/Logger.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;
using PLATO::Analysis::Encoding;
using PLATO::Utility::Logger;

using std::vector;
using std::string;
using std::fabs;
using std::log;
using std::exp;
using std::set;

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

	// We want to check for colinear columns in our data.  We will do that with
	// the SVD.  If no colinear columns are found, we couls always use that
	// to run our regression if we so choose (to save time)

	// NOTE: I don't want X to be const! hence the leading _
	_gsl_matrix_const_view X = gsl_matrix_const_view_array_with_tda(data, n_rows, n_cols, offset + n_cols);
	gsl_matrix* A = gsl_matrix_calloc(n_rows, n_cols);
	gsl_matrix* V = gsl_matrix_alloc(n_cols, n_cols);
	gsl_vector* S = gsl_vector_calloc(n_cols);
	gsl_vector* __ws_v = gsl_vector_alloc(n_cols);
	gsl_matrix* __ws_m = gsl_matrix_alloc(n_cols, n_cols);

	// P is our permutation matrix
	gsl_matrix* P = gsl_matrix_alloc(n_cols, n_cols);
	gsl_matrix_set_identity(P);

	gsl_matrix_memcpy(A, &X.matrix);

	gsl_linalg_SV_decomp_mod (A, __ws_m, V, S, __ws_v);

	// OK, now we look for all singular values less than machine epsilon
	// (we'll use float epsilon for double precision - lots of willge room there)
	unsigned int n_indep = n_cols;
	vector<unsigned int> idx_permu;
	while(gsl_vector_get(S,n_indep - 1) < std::numeric_limits<float>::epsilon()){
		// If we're here, the "n_indep - 1" column of V is a linear combination
		// of columns that adds to 0, so we want the last non-zero coefficient

		unsigned int p_idx = n_cols;

		// At the end of this (empty) loop, p_idx will be the index of the column
		// to drop
		while(--p_idx > 0 && fabs(gsl_matrix_get(V,p_idx,n_indep-1)) < std::numeric_limits<float>::epsilon());

		idx_permu.push_back(p_idx);
		--n_indep;
	}

	if(n_indep != n_cols){
		if(offset == 0){
			Logger::log_err("WARNING: colinear columns detected in linear regression, dropping "
							+ boost::lexical_cast<string>(n_cols - n_indep) + " parameters");
		}

		r->n_dropped = n_cols - n_indep;

		// OK, let's construct a permutation matrix by walking from the front of
		// the idx_permu vector and creating a temporary permutation matrix
		// and left-multiplying by P (which is now the identity)
		gsl_matrix* P_tmp = gsl_matrix_alloc(n_cols, n_cols);
		gsl_matrix* _P_work = gsl_matrix_calloc(n_cols, n_cols);
		for(unsigned int i=0; i<idx_permu.size(); i++){
			unsigned int p_idx = idx_permu[i];
			gsl_matrix_set_identity(P_tmp);

			// Now, we want to transpose the (p_idx) and (n_cols - i - 1) columns
			// so we set (p_idx, p_idx) = (n_cols-i-1,n_cols-i-1) = 0
			// and (p_idx,n_cols-i-1) = (n_cols-i-1,p_idx) = 1
			gsl_matrix_set(P_tmp, p_idx,p_idx, 0);
			gsl_matrix_set(P_tmp, n_cols-i-1,n_cols-i-1, 0);
			gsl_matrix_set(P_tmp, n_cols-i-1,p_idx, 1);
			gsl_matrix_set(P_tmp, p_idx,n_cols-i-1, 1);

			// Now, set P = P * P_tmp
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P, P_tmp, 0.0, _P_work);
			std::swap(P, _P_work);
		}
		gsl_matrix_free(P_tmp);
		gsl_matrix_free(_P_work);

		// Now, we're done with all of the SVD nonsense, we can set
		// A = X * P (remember, X is our original data, but it's const!)
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &X.matrix, P, 0.0, A);

		// Now, let's make X a view of only the first linearly independent columns
		// of A
		X = gsl_matrix_const_submatrix(A,0,0,n_rows, n_indep);
	}

	// At this point, we have our X matrix set appropriately
	// Let's set up our other data

	gsl_vector_const_view y_vec = gsl_vector_const_view_array(Y, n_rows);

	// First, let's make sure to initialize everything to 0, please!
	gsl_vector_view bv = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_vector_set_zero(&bv.vector);

	// now, get the proper sized array
	bv = gsl_vector_view_array(r->beta_vec, n_indep);
	gsl_vector* resid = gsl_vector_alloc(n_rows);

	// Again, let's make that covariance matrix all 0's!
	gsl_matrix* cov_mat = gsl_matrix_calloc(n_cols, n_cols);

	// We only want the leading n_indep x n_indep matrix here
	gsl_matrix_view cov_view = gsl_matrix_submatrix(cov_mat,0,0,n_indep, n_indep);
	gsl_matrix* cov = &cov_view.matrix;

	double chisq;

	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(n_rows, n_indep);

	// TODO: re-use the SVD we calculated previously if n_cols == n_indep
	gsl_multifit_linear(&X.matrix, &y_vec.vector, &bv.vector, cov, &chisq, ws);
	gsl_multifit_linear_residuals(&X.matrix, &y_vec.vector, &bv.vector, resid);

	// OK, now we have our coefficients!  It's that easy!
	// Well, sort of... now we have to un-permute everything if we had some issues
	// with colinearity

	// Note: to unpermute, multiply by P transpose!
	// Also, we need to unpermute both the rows AND columns of cov_mat
	bv = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_matrix* _cov_work = gsl_matrix_calloc(n_cols, n_cols);
	// permute columns
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, cov_mat, P, 0.0, _cov_work);
	// permute rows
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, _cov_work, 0.0, cov_mat);
	gsl_matrix_free(_cov_work);

	gsl_vector* _bv_work = gsl_vector_alloc(n_cols);
	gsl_vector_memcpy(_bv_work, &bv.vector);
	gsl_blas_dgemv(CblasTrans, 1.0, P, _bv_work, 0.0, &bv.vector);
	gsl_vector_free(_bv_work);

	r->coeffs.clear();
	r->stderr.clear();
	r->p_vals.clear();
	// add all the non-covariate coefficients

	// create a set of all of the removed indices
	set<unsigned int> permu_idx_set(idx_permu.begin(), idx_permu.end());

	for(unsigned int i=1+covar_names.size(); i<n_cols; i++){
		if(permu_idx_set.find(i) == permu_idx_set.end()){
			double c = gsl_vector_get(&bv.vector, i);
			double se = sqrt(gsl_matrix_get(cov_mat, i, i));
			r->coeffs.push_back(c);
			r->stderr.push_back(se);
			// t-val = | beta / stderr |
			if(encoding == Encoding::WEIGHTED){
				// this assumes that as df -> /inf, T -> Norm, and Norm^2 = ChiSq
				r->p_vals.push_back( gsl_cdf_chisq_Q( pow( c/se , 2) , 2));
			}else{
				r->p_vals.push_back(2*(1-gsl_cdf_tdist_P(fabs(c / se),n_rows-n_cols+1)));
			}
		} else {
			// If this is true, this column was dropped from analysis!
			r->coeffs.push_back(std::numeric_limits<float>::quiet_NaN());
			r->stderr.push_back(std::numeric_limits<float>::quiet_NaN());
			r->p_vals.push_back(std::numeric_limits<float>::quiet_NaN());
		}

	}

	addResult(r, null_result);


	double tss = gsl_stats_tss(Y, 1, n_rows);
	r->r_squared = 1 - chisq/tss;
	double null_rss = tss;
	unsigned int df = n_indep - reduced_vars;
	if (null_result){
		null_rss *= 1-null_result->r_squared;
		df += null_result->n_dropped;
	}


	double F=((null_rss - chisq) * (n_rows - n_indep))/(chisq * (df));

	// We want to see if there are extra degrees of freedom, which can happen in
	// the case of the "categorical" model
	// we have an extra df per marker in the categorically encoded model
	// We only have markers if we are not excluding markers and the
	// number of columns is at least as many as the number of covariates
	// (i.e. this isn;t the "null" model)
	unsigned int extra_df = (encoding == Encoding::WEIGHTED)
			* (!interactions || offset != 0)
			* (!exclude_markers) * (n_cols > covar_names.size() + 1)
			* (1 + pairwise);


	r->p_val = 1-gsl_cdf_fdist_P(std::max(0.0, F),df+extra_df,n_rows-n_indep-extra_df);

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

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(__ws_v);
	gsl_matrix_free(__ws_m);
	gsl_vector_free(resid);
	gsl_matrix_free(cov_mat);
	gsl_multifit_linear_free(ws);

	return r;
}

void LinearRegression::process(DataSet& ds){
	runRegression(ds);
}


}
}
