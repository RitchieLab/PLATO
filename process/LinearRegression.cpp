#include "LinearRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <set>
#include <algorithm>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>

#include "util/Logger.h"
#include "util/GSLUtils.h"

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
using std::pair;
using std::map;

namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string LinearRegression::stepname = LinearRegression::doRegister("linear");
const std::string LinearRegression::MPIName = LinearRegression::registerMPI("linear");


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

bool LinearRegression::initData(const PLATO::Data::DataSet& ds){

	vector<float>::const_iterator first_nonmiss = _pheno.begin();
	while(first_nonmiss != _pheno.end() && !std::isfinite(*first_nonmiss)){
		++first_nonmiss;
	}
	bool good_pheno = (first_nonmiss != _pheno.end());
	if(!good_pheno){
		Logger::log_err("ERROR: Given phenotype is completely missing", outcome_names.size() <= 1);
	}
	return good_pheno;
}

Regression::calc_fn& LinearRegression::getCalcFn() const{
	return (LinearRegression::calculate);
}

Regression::Result* LinearRegression::calculate(
		const double* Y, const double* data,
		unsigned int n_cols, unsigned int n_rows, unsigned int offset,
		unsigned int n_covars, const Regression::ExtraData* extra_data){

	//const extraData* extra_data = (const extraData*) (other_data);

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
		r->submodel = calculate(Y, data, reduced_vars, n_rows, offset + n_cols - (reduced_vars), new_covars, extra_data);

	}

	// We want to calculate the best fit for X*b = y

	gsl_matrix_const_view data_mv = gsl_matrix_const_view_array_with_tda(data, n_rows, n_cols, offset + n_cols);

	gsl_matrix* P = gsl_matrix_alloc(n_cols, n_cols);
	r->n_dropped = Utility::GSLUtils::checkColinear(&data_mv.matrix, P);

	gsl_matrix* A = gsl_matrix_alloc(n_rows, n_cols);

	// Now, we're done with all of the SVD nonsense, we can set
	// A = X * P (remember, X is our original data, but it's const!)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &data_mv.matrix, P, 0.0, A);

	unsigned int n_indep = n_cols - r->n_dropped;

	// Now, let's make X a view of only the first linearly independent columns
	// of A
	gsl_matrix_const_view X = gsl_matrix_const_submatrix(A,0,0,n_rows, n_indep);

	// At this point, we have our X matrix set appropriately
	// Let's set up our other data
	gsl_vector_const_view y_vec = gsl_vector_const_view_array(Y, n_rows);

	// First, let's make sure to initialize everything to 0, please!
	gsl_vector* beta = gsl_vector_calloc(n_cols);

	// now, get the proper sized array
	gsl_vector_view bv = gsl_vector_subvector(beta, 0, n_indep);
	gsl_vector* resid = gsl_vector_alloc(n_rows);

	// Again, let's make that covariance matrix all 0's!
	gsl_matrix* cov_mat = gsl_matrix_calloc(n_cols, n_cols);

	// We only want the leading n_indep x n_indep matrix here
	gsl_matrix_view cov_view = gsl_matrix_submatrix(cov_mat,0,0,n_indep, n_indep);
	gsl_matrix* cov = &cov_view.matrix;

	double chisq;

	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(n_rows, n_indep);

	gsl_multifit_linear(&X.matrix, &y_vec.vector, &bv.vector, cov, &chisq, ws);
	gsl_multifit_linear_residuals(&X.matrix, &y_vec.vector, &bv.vector, resid);

	// OK, now we have our coefficients!  It's that easy!
	// Well, sort of... now we have to un-permute everything if we had some issues
	// with colinearity

	// Note: to unpermute, multiply by P transpose!
	// Also, we need to unpermute both the rows AND columns of cov_mat
	//bv = gsl_vector_view_array(r->beta_vec, n_cols);
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
/*	gsl_vector* idx_permu = gsl_vector_calloc(n_cols);
	for(unsigned int i=0; i<n_cols; i++){
		gsl_vector_set(_bv_work, i, i);
	}
	gsl_blas_dgemv(CblasNoTrans, 1.0, P, _bv_work, 0.0, idx_permu);
	set<unsigned int> permu_idx_set(idx_permu->data + n_indep, idx_permu->data + n_cols);

	gsl_vector_free(idx_permu);
*/	gsl_vector_free(_bv_work);

	unsigned int idx_offset = 1+extra_data->base_covars;
	for(unsigned int i=0; i<n_cols-idx_offset; i++){
		double c = gsl_vector_get(beta, i+idx_offset);
		double se = sqrt(gsl_matrix_get(cov_mat, i+idx_offset, i+idx_offset));

		if(se > 0){

			r->coeffs[i] = c;
			r->stderr[i] = se;
			// t-val = | beta / stderr |
			if(extra_data->encoding == Encoding::WEIGHTED){
				// this assumes that as df -> /inf, T -> Norm, and Norm^2 = ChiSq
				r->p_vals[i] = gsl_cdf_chisq_Q( pow( c/se , 2) , 2);
			}else{
				r->p_vals[i] = 2*gsl_cdf_tdist_Q(fabs(c / se),n_rows-n_cols+1);
			}
		} else {
			// If this is true, this column was dropped from analysis!
			r->coeffs[i] = std::numeric_limits<float>::quiet_NaN();
			r->stderr[i] = std::numeric_limits<float>::quiet_NaN();
			r->p_vals[i] = std::numeric_limits<float>::quiet_NaN();
		}

	}

	// We want to see if there are extra degrees of freedom, which can happen in
	// the case of the "categorical" model
	// we have an extra df per marker in the categorically encoded model
	// We only have markers if we are not excluding markers and the
	// number of columns is at least as many as the number of covariates
	// (i.e. this isn;t the "null" model)

	// first, start off with just the number of independent columns
	r->df = findDF(P, reduced_vars, r->n_dropped, extra_data->extra_df_col, extra_data->n_extra_df_col);

	Result* curr_res = r;

	double tss = gsl_stats_tss(Y, 1, n_rows);
	string extraSuff = "";

	unsigned int df = 0;
	while(curr_res){
		df += curr_res->df;
		pair<float, float> pv_rsq = calcPVal(r, curr_res, chisq, tss, n_rows, df);

		if(curr_res == r){
			r->p_val = pv_rsq.first;
			r->r_squared = pv_rsq.second;
		} else if(curr_res->n_vars > 0) {
			extraSuff = boost::lexical_cast<string>(pv_rsq.first) + extra_data->sep;
			break;
		}

		curr_res = curr_res->submodel;
	}

	r->suffix += extraSuff;

	// I have no idea if the log_likelihood is correct!!
	r->log_likelihood = 0.5 * (-n_rows * (log(2*M_PI)+1 - log(n_rows) + log(chisq)));
	//double pv_test = gsl_cdf_chisq_Q(r->log_likelihood,1);

	// Make sure to clean up after yourself!

	gsl_vector_free(beta);
	gsl_matrix_free(A);
	gsl_matrix_free(P);
	gsl_vector_free(resid);
	gsl_matrix_free(cov_mat);
	gsl_multifit_linear_free(ws);

	return r;
}

pair<float, float> LinearRegression::calcPVal(Result* r, Result* curr_res, double chisq, double tss, unsigned int n_rows, unsigned int df){
	pair<float, float> pv_rsq;
	pv_rsq.first = 1;
	pv_rsq.second = 1 - chisq/tss;

	double null_rss = tss;
	if(curr_res->submodel){
		null_rss *= 1 - curr_res->submodel->r_squared;
	}

	if(df != 0){
		double F=((null_rss - chisq) * (n_rows - df))/(chisq * (df));
		pv_rsq.first = gsl_cdf_fdist_Q(std::max(0.0, F),df,n_rows-df);
	}

	return pv_rsq;

}

void LinearRegression::process(DataSet& ds){
	runRegression(ds);
}

pair<unsigned int, const char*> LinearRegression::calculate_MPI(unsigned int bufsz, const char* buf){
	return Regression::calculate_MPI(bufsz, buf, LinearRegression::calculate);
}

}
}
