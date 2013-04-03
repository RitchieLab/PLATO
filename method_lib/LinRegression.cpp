//LinRegression.cpp
#include "LinRegression.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
#include <cmath>;

namespace Methods{

///
/// Constructor
///
LinRegression::LinRegression(){
  set=NULL;
  initialize();
}


///
/// Alternative constructor
///
LinRegression::LinRegression(DataSet* ds){
  initialize();
  resetDataSet(ds);
}


///
/// Initialize starting variables
///
void LinRegression::initialize(){
  snp_interaction  = 0;
  r2 = adjusted_r2 = likelihood = lrt_pval = 0.0;
  f_pval = 1.0;

  ModelTypeMap["DOMINANT"] = Dominant;
  ModelTypeMap["RECESSIVE"] = Recessive;
  ModelTypeMap["ADDITIVE"] = Additive;
  
  modType = Additive;
  missingValue = 3;
  missingCoValue = -99999;
  
  set_model();
  
}

///
/// run linear regression on the indicated loci
/// @param loci 
/// @param set DataSet to use for this calculation
///
void LinRegression::calculate(vector<unsigned int>& loci, DataSet* set){
  resetDataSet(set);
  calculate(loci);
}

///
/// runs linear regression on the dataset set within the method
/// @param loci 
///
void LinRegression::calculate(vector<unsigned int>& loci){
  vector<unsigned int> covars;
  calculate(loci, covars);
}

///
/// runs linear regression on combination of traits, loci, and covariates
/// @param loci
/// @param covars Vector containing indices for covariates to include in model 
///
void LinRegression::calculate(vector<unsigned int>& loci, vector<unsigned int>& covars){
  prepare_input(loci, covars);     
}


///
/// Set up GSL structures for running linear regression.
/// @param loci vector pass empty if no markers included in analysis
/// @param covars vector pass empty if no covariates included in analysis
///
void LinRegression::prepare_input(vector<unsigned int>& loci, vector<unsigned int>& covars){

  unsigned int numInds = set->num_inds();
  
  // convert loci indexes to marker map
  vector<unsigned int> converted_loci = convert_loc_order(loci);
  
  // set reference alleles (needed for different models)
  vector<int> ref_alleles(converted_loci.size(), 0);
  for(unsigned int curr_mark=0; curr_mark < converted_loci.size(); curr_mark++)
  {
    ref_alleles.at(curr_mark) = set->get_locus(loci.at(curr_mark))->getReferentIndex();
  }  
  
  unsigned int num_covars = covars.size();
  unsigned int num_loci = loci.size();
  
  int n_vars =  num_loci + num_covars + snp_interaction;
  
  vector<vector<double> > analysis_matrix;
  
  unsigned int currInd;
  // first column will be for the outcome variable
  vector<double> row(n_vars+1,0.0);
  int curr_column;
  bool any_missing;

  for(currInd=0; currInd < numInds; currInd++)
  {
    curr_column=0;
    any_missing=false;
    // ignore individuals that are not enabled  
    if(!(*set)[currInd]->isEnabled())
		  continue;
		
		// add phenotype from indicated trait
		row.at(curr_column++) = Y[currInd];

    for(unsigned int i=0; i < num_loci; i++)
    {
      if((*set)[currInd]->get_genotype(converted_loci.at(i)) != missingValue)
        row.at(curr_column++) = geno_convert.at(ref_alleles.at(i)).at((*set)[currInd]->get_genotype(converted_loci.at(i)));
// row.at(curr_column++)=(*set)[currInd]->get_genotype(converted_loci.at(i));
      else
      {
        any_missing = true;
        break;
      }
    }
    
    // add interaction when needed 
    if(snp_interaction){
      row.at(curr_column) = row[curr_column-1] * row[curr_column-2];
      curr_column++;
    }
		
	  for(unsigned int i=0; i<num_covars; i++){
	    if((*set)[currInd]->getCovariate(covars.at(i)) != missingCoValue)
        row.at(curr_column++) = (*set)[currInd]->getCovariate(covars.at(i));
      else{
        any_missing = true;
        break;
      } 
	  }

    
    if(!any_missing)
      analysis_matrix.push_back(row);
  }
  
  // total number of cases will equal the number of analysis matrix rows
  ngenotypes = int(analysis_matrix.size());
  
  /// now that the analysis matrix has been constructed can do calculation
  calculate_linreg(analysis_matrix);
}


///
/// Utilizes GSL library to run core analysis and then fills other parameters
/// for return as needed by user.
/// @param analysis_matrix Each row corresponds to an individual in the set. The first
/// column should be the outcome variable.
///
void LinRegression::calculate_linreg(vector<vector<double> >& analysis_matrix){
  gsl_matrix *X, *cov;
  gsl_vector *y, *c;
  double chisq;
  
  int n = int(analysis_matrix.size());
  // number of variables is number - the outcome variable but with an additional 
  // constant of 1.0 for the intercept in the first position
  int n_vars = int(analysis_matrix[0].size()); 
  
  X = gsl_matrix_alloc(n, n_vars);
  y = gsl_vector_alloc(n);
  c = gsl_vector_alloc(n_vars);
  cov = gsl_matrix_alloc(n_vars, n_vars);
  
  double* y_array = new double[n]; // for TSS calculation
  
  for(int i=0; i<n; i++){
    gsl_vector_set(y,i,analysis_matrix[i][0]); // set outcome
    gsl_matrix_set(X,i,0,1.0); // intersection term
    y_array[i] = analysis_matrix[i][0];
    for(int j=1; j < n_vars; j++){
      gsl_matrix_set(X,i,j,analysis_matrix[i][j]);
    }
  }
  
  // create GSL workspace
  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(n,n_vars);
  // run analysis
  gsl_multifit_linear(X,y,c,cov,&chisq, work);
  // free GSL workspace
  gsl_multifit_linear_free (work);
  
  // calculate and fill results
  coefficients.clear();
  std_errors.clear();
  tt_vals.clear();
  coeff_pvals.clear();
  
  #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
  
  int df = n-n_vars; // degrees of freedom
  
  double stderr, tval, pval;

  coeff_intercept = gsl_vector_get(c,0);  
  
  for(int i=1; i<n_vars; i++){
    coefficients.push_back(gsl_vector_get(c,i));
    stderr = sqrt(COV(i,i));
    tval = gsl_vector_get(c,i)/stderr;
//     pval = tval<0?2*(1-gsl_cdf_tdist_P(-tval,n-4)):2*(1-gsl_cdf_tdist_P(tval,n-4));
    pval = tval<0?2*(1-gsl_cdf_tdist_P(-tval,df)):2*(1-gsl_cdf_tdist_P(tval,df));
    
    std_errors.push_back(stderr);
    tt_vals.push_back(tval);
    coeff_pvals.push_back(pval);
  }
  
  // calculate r-squared and adjusted r-squared
  double tss = gsl_stats_tss(y_array,1,n);
  
  r2 = 1-(chisq/tss);
  adjusted_r2 = 1-double(n-1)/df*(1-r2);
  
  // calculate F statistic
  double ymean = gsl_stats_mean(y_array,1,n);
  
  #define VAL(i,j) (gsl_matrix_get(X,(i),(j)))
  // sum difference between predicted Y and mean
  // sum difference between predicted Y and actual Y
  double ymeansum=0.0;
  double yactualsum=0.0;
  vector<double> residuals;
  double predicted_y;
  for(int i=0; i<n; i++){
    predicted_y = gsl_vector_get(c,0);
    for(int j=1; j<n_vars; j++){
      predicted_y += gsl_vector_get(c,j)*VAL(i,j);
    }
    residuals.push_back(gsl_vector_get(y,i) - predicted_y);
    double ymeandiff = predicted_y - ymean;
    ymeansum += ymeandiff * ymeandiff;
    double yactualdiff = predicted_y - y_array[i];
    yactualsum += yactualdiff * yactualdiff;
  }
  
  double p = n_vars-1;
  double fnumerator = ymeansum /p;
  double fdenominator = yactualsum / df;
  double fstatistic = fnumerator / fdenominator;
  
  // fstat p-value can be used as measure of overall model p-value
  f_pval = 1-gsl_cdf_fdist_P(fstatistic,p,df);
  
  // calculate the log-likelihood
  // calculate squared sum of residuals
  vector<double>::iterator resiter;
  double sumressq=0.0;
  for(resiter = residuals.begin(); resiter != residuals.end(); resiter++){
    sumressq += *resiter * *resiter;
  }  
  double log_sum_ressq = log(sumressq);
  double N = residuals.size();
  likelihood = 0.5 * (-N * (log(2*M_PI)+1 - log(N) + log_sum_ressq)); 
  
  lrt_pval = 1-gsl_cdf_chisq_P(likelihood,1);
  
  delete [] y_array;
  
  // free all GSL vectors and matrices
  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (c);
  gsl_matrix_free (cov);  
  
}




///
/// Returns encoding for a genotype based on the model set.
/// @param geno Genotype value being checked
/// @param referent_allele Index of the allele that is the referent for the marker
/// @return recoded value for the genotype
///
int LinRegression::get_geno_conversion(int geno, int referent_allele){
	return geno_convert.at(referent_allele).at(geno);
}



///
/// Indexes refer to the order of the markers in map location.
/// To use in sample need to convert them to the indexes in the Samples.
/// @param loci vector of unsigned int that will be converted
/// @return returns vector with updated indexes
///
vector<unsigned int> LinRegression::convert_loc_order(vector<unsigned int>& loci){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  vector<unsigned int>::iterator iter;
  vector<unsigned int> converted_indexes;

  for(iter=loci.begin(); iter!=loci.end(); iter++){
    converted_indexes.push_back((*markers)[*iter]->getLoc());
  }

  return converted_indexes;
}


///
/// Sets the parameters using StepOptions class
/// @param options StepOptions containing options
///
void LinRegression::set_parameters(StepOptions* o){
  options = o;
  setModelType(options->getLRModelType());
  
  
  setDependent();
  
}


//
/// Sets the model type to use in the calculation
/// @param modelType string containing model type to run
/// @throws Exception on error
///
void LinRegression::setModelType(string modelType){
  if(ModelTypeMap.find(modelType) != ModelTypeMap.end()){
    modType = ModelTypeMap[modelType];
  }
  else{
    throw MethodException(modelType + " is not a valid LR model type");
  }
  set_model();
}

///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void LinRegression::resetDataSet(DataSet* ds){
  set = ds;
  missingValue = set->get_missing_value();
  missingCoValue = set->get_missing_covalue();
  markers = set->get_markers();
}


///
/// Return trait index to use if set
///
///
void LinRegression::setDependent() {

  Y.clear();
  if(options->getUsePheno()){
    int index = options->getPhenoLoc();
    if (options->getPhenoName() != "") {
		  index = set->get_trait_index(options->getPhenoName());
		}   
		for (int i = 0; i < set->num_inds(); i++){
		  Y.push_back(set->get_sample(i)->getPheno(index));
		}
  }
  else{
    for (int i = 0; i < set->num_inds(); i++){
      Y.push_back(set->get_sample(i)->getPheno());
    }
  }
  
}


///
/// Sets values to use for each genotype based on model selected <br>
/// Only works for SNPs currently.
///
void LinRegression::set_model(){
  // set up 2-D array
  geno_convert.assign(2, vector<unsigned int>(4,0));

  switch(modType){
    case Dominant:
      geno_convert.at(0).at(0) = 1;
      geno_convert.at(0).at(1) = 1;
      geno_convert.at(0).at(3) = 3;
      geno_convert.at(1).at(1) = 1;
      geno_convert.at(1).at(2) = 1;
      geno_convert.at(1).at(3) = 3;
      break;
    case Recessive:
      geno_convert.at(0).at(0) = 1;
      geno_convert.at(0).at(3) = 3;
      geno_convert.at(1).at(2) = 1;
      geno_convert.at(1).at(3) = 3;
      break;
    case Additive:
      geno_convert.at(0).at(1) = 1;
      geno_convert.at(0).at(0) = 2;
      geno_convert.at(0).at(3) = 3;
      geno_convert.at(1).at(1) = 1;
      geno_convert.at(1).at(2) = 2;
      geno_convert.at(1).at(3) = 3;
      break;
  }
}

} // end namespace 