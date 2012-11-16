//CondLogRegression.h

#ifndef __COND_LOGISTICREGRESSION_H__
#define __COND_LOGISTICREGRESSION_H__
#include "DataSet.h"
#include "FlatIndex.h"
#include "StepOptions.h"
#include "conditional_lr.h"
#include "MethodException.h"

namespace Methods{
///
/// Currently only works on datasets that 
/// have only one affected and one unaffected
/// individual per pedigree
///

/// Method class that performs conditional logistic regression
class CondLogRegression{

  public:
    CondLogRegression();
    CondLogRegression(DataSet* ds);
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);

    /// sets type of model to use in LR calculation
    void setModelType(string modType);
    
    /// run logistic regression on the indicated loci
    void calculate(vector<unsigned int>& loci, DataSet* set);
  
    /// runs logistic regression on the dataset set within the method
    void calculate(vector<unsigned int>& loci);
  
    /// returns the overall model p value
    double getOverallP(){return p_value;}
    
    /// returns the log likelihood ratio
    double getLogLikelihoodRatio(){return likelihood_ratio;}
    
    /// returns the coefficients for the terms
    vector<double> getCoefficients(){return coefficients;}
    
    /// returns the covariates
    vector<double> getCovariates(){return covariates;}
    
    /// specifies model type for use in this run of logistic regression
    enum ModelTypes{
      Dominant,
      Recessive,
      Additive
    };

    /// sets whether to include full interaction term in analysis (false excludes full interaction)
    void setFullInteraction(bool interact){fullInteraction = interact;}

    /// if true, interaction terms are included in a multilocus model
    void setIncludeInteractions(bool interact){interaction_included = interact;}
  
    /// sets maximum number of iterations that algorithm will use in trying to reach convergence
    void setMaximumIterations(unsigned int maxIter){maxIterations = maxIter;}
    /// returns current number of iterations
    unsigned int getMaximumIterations(){return maxIterations;}
    
  private:
  
    DataSet* set;
  
    struct conditional_lr_parameters{
      doublereal * z;
      integer * nca;
      integer * nct;
      integer * ivar;
      integer * ins;
      doublereal * b;
      integer ns, nimax, nmax, nvmax, nvmax1, nv, nmax1;
      doublereal *cntr, *w, *wb, *wdb, *wd2b, *u,  *db, *d2b, *dl, *cov, *covi;
      doublereal chi2, st;
      integer ifault;      
    };  
  
    void initialize();
    void initialize_interactions();
    void set_model();
    void select_individuals();
    float run_conditional_lr(conditional_lr_parameters& param);
    void allocate_arrays(conditional_lr_parameters& params, vector<unsigned int>& loci);
    void adjust_arrays(conditional_lr_parameters& params, vector<unsigned int>& loci);
    void fill_arrays(conditional_lr_parameters& params);
    void fill_params(conditional_lr_parameters& params, vector<unsigned int>& loci);
    void output_params(conditional_lr_parameters& params);
    void fill_z_array(conditional_lr_parameters& params, vector<unsigned int>& loci);
    void free_params(conditional_lr_parameters& params);
    double ChiSq(double x, unsigned int n);
    double norm(double z);
    
    float likelihood_ratio, p_value;
    vector<double> covariates, coefficients;
    conditional_lr_parameters param;

    vector<vector<Sample*> > strata_inds;
    vector<unsigned int> include_strata;

    bool fullInteraction, interaction_included;
  
    ModelTypes modType;
    map<string, ModelTypes> ModelTypeMap;
    vector<vector<vector<unsigned int> > > interaction_lists;

    unsigned int maxLocusValue, missingValue, maxIterations, modelSize, 
      LociComboMin, LociComboLimit, defaultComboInterval, original_strat_size;
    
//    vector<unsigned int> geno_convert;
    vector<vector<unsigned int> >geno_convert;

    double PiD2;    

};
};

#endif
