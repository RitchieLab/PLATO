//LinRegression.h

#ifndef __LINREGRESSION_H__
#define __LINREGRESSION_H__

// #include "DataSet.h"
#include "Regression.h"
// #include "StepOptions.h"
#include "MethodException.h"
#include <set>

///
/// Linear regression in method library. Utilizes GSL
/// to perform main analysis and then calculates
/// remaining values for use in plato.
///

/// Method class that performs linear regression
namespace Methods{
class LinRegression: public Regression{

  public:
    LinRegression();
    LinRegression(DataSet* ds);
  
    /// run linear regression on the indicated loci
    void calculate(vector<unsigned int>& loci, DataSet* set);

    /// runs linear regression on the dataset set within the method
    void calculate(vector<unsigned int>& loci);

    /// runs linear regression on combination of traits, loci, and covariates
    void calculate(vector<unsigned int>& loci, vector<unsigned int>& covars);

    /// specifies model type for use in this run of linear regression
    enum ModelTypes{
      Dominant,
      Recessive,
      Additive
    };

    /// sets type of model to use in LR calculation
    void setModelType(string modType);
    
    /// sets DataSet
    void resetDataSet(DataSet* ds);

    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);

    /// Get model coefficients
    vector<double> getCoefficients(){return coefficients;}
    
    /// Get standard errors of coefficients
    vector<double> getCoeffStandardErr(){return std_errors;}
    
    /// Get p values for coefficients
    vector<double> getCoeffPValues(){return coeff_pvals;}

    /// returns intersect coefficient
    double getCoeffIntercept(){return coeff_intercept;}

    /// Get p value for model
    double getOverallP(){return f_pval;}
    
    /// Get overall score
    double getOverallScore(){return f_statistic;}
    
    /// Get r-squared for model
    double getRSquared(){return r2;}
    
    /// Get adjusted r-squared for model
    double adjusted_rsquared(){return adjusted_r2;}
    
    /// Get log-likelihood for model
    double getLLR(){return likelihood;}
    
    /// return # of genotypes used in calculation
    int getNumGenotypes(){ return ngenotypes;}
    
    /// sets whether to include interaction term for SNPs
    void setIncludeInteractions(bool include){
      if(include)
        model_interaction=1;
      else
        model_interaction=0;
    }
    
  private:

    void prepare_input(vector<unsigned int>& loci, vector<unsigned int>& covars);
    void calculate_linreg(vector<vector<double> >& analysis_matrix);
    vector<unsigned int> convert_loc_order(vector<unsigned int>& loci);
    int get_geno_conversion(int geno, int referent_allele);
    void set_model();
    void initialize();
    void setDependent();
  
    StepOptions * options;
    DataSet* set;
    int n_vars, n_inds, model_interaction, ngenotypes;
    vector<double> coefficients, std_errors, tt_vals, coeff_pvals;
    double f_pval, r2, adjusted_r2, likelihood, lrt_pval, missingCoValue, coeff_intercept,
    	f_statistic;
    unsigned int missingValue;
    
    ModelTypes modType;
    map<string, ModelTypes> ModelTypeMap;
    vector<vector<unsigned int> >geno_convert;
    vector<Marker*> * markers;
    std::set<int> skipInd;
    vector<double> Y;

};
};
#endif