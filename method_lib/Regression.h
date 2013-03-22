//Regression.h

#ifndef __REGRESSION_H__
#define __REGRESSION_H__

#include "DataSet.h"
#include "StepOptions.h"

///
/// Abstract base class for regression methods
///

/// Method class that performs logistic regression
namespace Methods{
class Regression{

  public:
    Regression(){};
    Regression(DataSet* ds){};

    /// run regression on the indicated loci using the DataSet passed
    virtual void calculate(vector<unsigned int>& loci, DataSet* set)=0;

    /// run regression on the indicated 
    virtual void calculate(vector<unsigned int>& loci)=0;

    /// returns the overall model p value
    virtual double getOverallP()=0;

    /// returns the log likelihood ratio
    virtual double getLLR()=0;

    /// returns the coefficient p value of the full interaction term
//     virtual double getFullInteractionP()=0;

    /// returns the coefficients for the terms
    virtual vector<double> getCoefficients()=0;

    /// returns the standard errors for the coefficients
    virtual vector<double> getCoeffStandardErr()=0;

    /// return # of genotypes used in calculation
    virtual int getNumGenotypes()=0;

    /// returns the p-values for coefficients
    virtual vector<double> getCoeffPValues()=0;

    /// returns intersect coefficient
    virtual double getCoeffIntercept()=0;

    /// include or don't include highest order interaction of model
//     virtual bool setIncludeInteraction(bool inter)=0;

    /// set parameters for method using StepOptions class
    virtual void set_parameters(StepOptions* options)=0;

    /// sets type of model to use in LR calculation
    virtual void setModelType(string modType)=0;

    /// returns adjusted rsquared (or pseudo-rsquared for logistic regression)
    virtual double adjusted_rsquared()=0;

    /// sets DataSet
    virtual void resetDataSet(DataSet* ds)=0;
    
    virtual void setIncludeInteractions(bool include)=0;
    
//     virtual double getRSquared()=0;


  private:

};
};
#endif
