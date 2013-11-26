//BoxCoxTransform.h

#ifndef __BOXCOXTRANSFORM_H__
#define __BOXCOXTRANSFORM_H__

#include "DataSet.h"
#include "StepOptions.h"
#include "MethodException.h"

///
/// Transforms continuous data using the Box Cox transformation.  Values
/// can be converted back to the original score by using the same object or storing
/// the lambda values and then providing them to a new instance of the class. Values
/// that match the designated missing value are ignored.
/// Based on pg 417-419 in Biometry by Robert R. Sokal/F. James Rohlf (Third Edition)
///

namespace Methods{

class BoxCoxTransform{

  public:
    /// Transforms all data in set.
    void TransformData(DataSet* set);
  
    /// Transforms a covariate
    void TransformCovar(DataSet* set, int covar);

    /// Transforms a trait
    void TransformTrait(DataSet* set, int trait);
  
    /// Returns all data back to original values
    void UndoTransform(DataSet* set);
          
    /// Sets lambda values for reversing transformation
    void SetLambdas(vector<double>& covariate_lambdas, vector<double>& trait_lambdas){
      covar_lambda=covariate_lambdas;
      trait_lambda=trait_lambdas;
    }
     
    /// Undoes transformation of a covariate
    void UndoCovariate(DataSet* set, int covar);
  
    /// Undoes transformation of a trait
    void UndoTrait(DataSet* set, int trait);
    
    /// Returns covariate lambda values 
    inline vector<double> GetCovarLambda(){return covar_lambda;}
    
    /// Returns trait lambda values
    inline vector<double> GetTraitLambda(){return trait_lambda;}
    
  private:
  
    double CalculateLikelihood(vector<float>& values, float log_total,
      double lambda);

    double OptimizeLambda(vector<float>& values, double log_total, 
      double startlam, double interval);

    double GetLambda(vector<float>& values);
  
    double IterateBest(vector<float>& values, double log_total, double left, double right);
  
    double Transform(double val, double lambda);
    double Undo(double val, double lambda);
    vector<double> covar_lambda, trait_lambda;
  
};

}


#endif
