//LogisticRegression.h

#ifndef __LOGISTICREGRESSION_H__
#define __LOGISTICREGRESSION_H__

#include "DataSet.h"
#include "FlatIndex.h"
#include "StepOptions.h"
#include "MethodException.h"

///
/// Calculates logistic regression in method library
/// used by plato and wasp.  This method currently
/// only accepts snps for markers.
///

/// Method class that performs logistic regression
namespace Methods{
class LogisticRegression{

  public:
    LogisticRegression();
    LogisticRegression(DataSet* ds);
  
    /// run logistic regression on the indicated loci
    void calculate(vector<unsigned int>& loci, DataSet* set);

    /// runs logistic regression on the dataset set within the method
    void calculate(vector<unsigned int>& loci);
    
    /// runs logistic regression on combination of traits, loci, and covariates
    void calculate(vector<unsigned int>& loci, vector<unsigned int>& covars, 
      vector<unsigned int> & traits);
  
    /// returns the overall model p value
    double getOverallP(){return overallPvalue;}
    
    /// returns the log likelihood ratio
    double getLLR(){return LLR;}
    
    /// returns the coefficient p value of the full interaction term
    double getFullInteractionP(){return coeffPvalue;}
    
    /// returns the coefficients for the terms
    vector<double> getCoefficients(){return coefficients;}

    /// returns the standard errors for the coefficients
    vector<double> getCoeffStandardErr(){return standard_errors;}
    
    /// returns intersect coefficient
    double getCoeffIntercept(){return coeff_intercept;}
  
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);
  
    /// specifies model type for use in this run of logistic regression
    enum ModelTypes{
      Dominant,
      Recessive,
      Additive
    };

    /// sets whether to include full interaction term in analysis (false excludes full interaction)
    void setFullInteraction(bool interact){fullInteraction = interact;}

    /// if true, interaction terms are included in a multilocus model
    void setIncludeInteractions(bool interact){includeInteractions = interact;}

    /// sets type of model to use in LR calculation
    void setModelType(string modType);
  
    /// sets maximum number of iterations that algorithm will use in trying to reach convergence
    void setMaximumIterations(unsigned int maxIter){maxIterations = maxIter;}
    /// returns current number of iterations
    unsigned int getMaximumIterations(){return maxIterations;}
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);

	  /// gets conversion value of genotype based on model and which allele is the referent allele
	  int get_geno_conversion(int geno, int ref_allele_index);
    
  private:
    DataSet* set;

    void initialize();
    void initialize_interactions();
    void set_model();
    double norm(double z);
    void summarize_data(vector<unsigned int> & genos);
    void initialize_summary(unsigned int currModelSize);
    unsigned int ix(int j,int k,int nCols);
    void zero_summary(unsigned int array_size, unsigned int model_size);
    double ChiSq(double x, unsigned int n);
    
    void calculateLR(vector<vector<double> >& data, bool summary_data, vector<int>& includedCells);
    vector<unsigned int> convert_loc_order(vector<unsigned int>& loci);
    
    vector<double> coefficients, standard_errors;
    double coeffPvalue, LLR, overallPvalue, coeff_intercept, missingCoValue;
    bool fullInteraction, includeInteractions;
  
    ModelTypes modType;
    map<string, ModelTypes> ModelTypeMap;
  
    FlatIndex indexConverter;
    vector<vector<int> > includedIndexes;
    unsigned int maxLocusValue, missingValue, maxIterations, modelSize, 
      LociComboMin, LociComboLimit, defaultComboInterval;
    
    vector<vector<unsigned int> >geno_convert;
    vector<Marker*> * markers;

    vector<vector<double> > summary_data;
    vector<vector<vector<unsigned int> > > interaction_lists;
    double PiD2;    

};
};
#endif
