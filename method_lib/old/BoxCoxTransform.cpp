//BoxCoxTransform.cpp

#include "BoxCoxTransform.h"

namespace Methods{

///
/// Transforms all data in set.
/// @param set DataSet
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::TransformData(DataSet* set){

  int nCovars = set->num_covariates();
  // reset the lambda vectors to hold the values
  covar_lambda.assign(nCovars,1.0);
  
  for(int covar=0; covar < nCovars; covar++){
    TransformCovar(set, covar);
  }
  
  int nTraits = set->num_traits();
  trait_lambda.assign(nTraits,1.0);
  for(int trait=0; trait < nTraits; trait++){
    TransformTrait(set, trait);
  }
}


///
/// Transforms a covariate.
/// @param set DataSet
/// @param covar Index of covariate to transform
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::TransformCovar(DataSet* set, int covar){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  
  // covariate lambda array hasn't been set assign it to be correct size
  if(int(covar_lambda.size()) != set->num_covariates()){
    // reset the lambda vectors to hold the values
    covar_lambda.assign(set->num_covariates(),1.0);  
  }
  
  // fill vector for calling GetLambda
  vector<float> covars;

  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getCovariate(covar) != missingCoValue){
      covars.push_back(samp->getCovariate(covar));
    }
  }

  double lambda = GetLambda(covars);

  covar_lambda.at(covar) = lambda;
  
  // set covariates
  int c=0;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getCovariate(covar)!=missingCoValue){
      samp->setCovariate(Transform(covars[c++], lambda), covar);
    }
  }

}


///
/// Calculates and returns lambda.  Method is based on the book
/// Applied Regression Analysis by Draper and Smith (3rd edition)
/// pages 280-281.
/// @param values Contains original values
/// @return lambda
///
double BoxCoxTransform::GetLambda(vector<float>& values){
  
  double log_total=0.0;
  for(unsigned int i=0; i<values.size(); i++){
    log_total += log(values.at(i));
  }

  double interval = 1.0;
  
  double start_lam = 0;
  
  // when find a point that is descending 
  return OptimizeLambda(values, log_total, start_lam, interval);
}



///
/// Performs hill climb to find the lambda
/// @param values
/// @param log_total
/// @param startlam
/// @param interval
/// @return lambda
///
double BoxCoxTransform::OptimizeLambda(vector<float>& values, double log_total, 
  double startlam, double interval){
  
  double lambda = startlam + interval;
  double last_lambda = startlam;
  double lastLL = CalculateLikelihood(values, log_total, startlam);
  double LL;
  
  while(fabs(last_lambda - lambda) > 0.0001){
    
    LL = CalculateLikelihood(values, log_total, lambda);

    if(LL < lastLL){
      // adjust interval to start moving in opposite direction
      interval = -interval/10;
    }
    
    lastLL = LL;
    last_lambda = lambda;
    lambda += interval;
  }
  
  return lambda;
  
}


///
/// Examines values from the left lambda to the right lambda and
/// finds the high point using a binary search
/// @param left Left lambda 
/// @param right Right lambda
/// @return best lambda value
///
double BoxCoxTransform::IterateBest(vector<float>& values, double log_total,
  double left, double right){
  
  double last_lambda = right;
  double lambda = left + (right - left) / 2;
  double LL, rightLL, leftLL;
  
  leftLL = CalculateLikelihood(values, log_total, left);
  rightLL = CalculateLikelihood(values, log_total, right);
  
  int nIterations = 0;

  while(fabs(last_lambda - lambda) > .001 && nIterations < 200){

    LL = CalculateLikelihood(values, log_total, lambda);
    if(LL > rightLL){ // set midpoint to be right side of search
      right = lambda;
      rightLL = LL;
    }
    else{ // set midpoint to be left side of search
      left = lambda;
      leftLL = LL;
    }
    
    last_lambda = lambda;
    lambda = left + (right - left) / 2;
    nIterations++;
  }

  return lambda;
}



///
/// Calculates likelihood for the data passed
/// @param vals 
/// @param log_total natural log total of the data
/// @param lambda
/// @return lambda
///
double BoxCoxTransform::CalculateLikelihood(vector<float>& values, float log_total,
  double lambda){
  
  int n = int(values.size());  // number of samples
  int v = n-1; // degrees of freedom

  vector<double> transformed(n, 0);
  
  // first term depends on the variance of transformed values (s2t)
  double s2t=0.0, transmean = 0.0;
  for(int i=0; i<n; i++){
    if(fabs(lambda) > 0.000001)
      transformed.at(i) = (pow(values.at(i), lambda)-1)/lambda;
    else
      transformed.at(i) = log(values.at(i));
    transmean += transformed.at(i);
  }
  transmean /= n;
  
  for(int i=0; i<n; i++){
    s2t += pow(transformed.at(i)-transmean, 2);
  }
  s2t = s2t / v;
  
  // calculate and return likelihood
  return -double(v)/2 * log(s2t) + (lambda-1) * (double(v)/n) * log_total;
  
}



///
/// Transforms a trait.
/// @param set DataSet
/// @param trait Index of trait to transform
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::TransformTrait(DataSet* set, int trait){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;

  if(int(trait_lambda.size()) != set->num_traits()){  
    trait_lambda.assign(set->num_traits(),1.0);
  }
  
  // fill vector for calling GetLambda
  vector<float> traits;

  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getTrait(trait) != missingCoValue){
      traits.push_back(samp->getTrait(trait));
    }
  }
  
  double lambda = GetLambda(traits);
  // covars contains the new transformed values
  trait_lambda.at(trait) = lambda;
  
  // set covariates
  int t=0;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getTrait(trait)!=missingCoValue)
      samp->setTrait(Transform(traits[t++], lambda), trait);
  }
}

  
/// 
/// Returns all data back to original values
/// @param set DataSet
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::UndoTransform(DataSet* set){
  int nCovars = set->num_covariates();
  for(int covar=0; covar < nCovars; covar++){
    UndoCovariate(set, covar);
  }
  
  int nTraits = set->num_traits();
  for(int trait=0; trait < nTraits; trait++){
    UndoTrait(set, trait);
  }  
}

///
/// Undoes transformation of a covariate
/// @param set DataSet
/// @param covar Index of covariate to undo
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::UndoCovariate(DataSet* set, int covar){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getCovariate(covar)==missingCoValue)
      continue;
    samp->setCovariate(Undo(samp->getCovariate(covar), covar_lambda.at(covar)), covar);
  }  
}



///
///  Undoes transformation of a trait
/// @param set DataSet
/// @param trait Index of trait to undo
/// @throws MethodException when negative value encountered
///
void BoxCoxTransform::UndoTrait(DataSet* set, int trait){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getTrait(trait)==missingCoValue)
      continue;
    samp->setTrait(Undo(samp->getTrait(trait), trait_lambda.at(trait)), trait);
  }
}


///
/// Transforms value and returns transformed value
/// @param val
/// @param lambda
/// @return transformed value
/// @throws MethodException when value is negative
///
double BoxCoxTransform::Transform(double val, double lambda){
  if(val < 0)
    throw MethodException("Negative value passed to BoxCox Transform");
  if(fabs(lambda) > 0.000001)
    return (pow(val, lambda)-1)/lambda;  
  else
    return log(val);
}

///
/// Undoes transformation and returns original value
/// @param val
/// @param lambda
/// @return original value
/// @throws MethodException when value is negative
///
double BoxCoxTransform::Undo(double val, double lambda){
  
  if(fabs(lambda) > 0.000001)
    return pow((val * lambda)+1, 1/lambda);
  else
    return exp(val);
  
}


}

