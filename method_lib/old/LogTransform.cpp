//LogTransform.cpp

#include "LogTransform.h"

namespace Methods{

///
/// Transforms all data in set.
/// @param set DataSet
/// @throws MethodException when negative value encountered
///
void LogTransform::TransformData(DataSet* set){
  int nCovars = set->num_covariates();
  for(int covar=0; covar < nCovars; covar++){
    TransformCovar(set, covar);
  }
  
  int nTraits = set->num_traits();
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
void LogTransform::TransformCovar(DataSet* set, int covar){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getCovariate(covar)==missingCoValue)
      continue;
    samp->setCovariate(Transform(samp->getCovariate(covar)), covar);
  }
}



///
/// Transforms a trait.
/// @param set DataSet
/// @param trait Index of trait to transform
/// @throws MethodException when negative value encountered
///
void LogTransform::TransformTrait(DataSet* set, int trait){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getTrait(trait)==missingCoValue)
      continue;
    samp->setTrait(Transform(samp->getTrait(trait)), trait);
  }
}

  
/// 
/// Returns all data back to original values
/// @param set DataSet
/// @throws MethodException when negative value encountered
///
void LogTransform::UndoTransform(DataSet* set){
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
void LogTransform::UndoCovariate(DataSet* set, int covar){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getCovariate(covar)==missingCoValue)
      continue;
    samp->setCovariate(Undo(samp->getCovariate(covar)), covar);
  }  
}


///
///  Undoes transformation of a trait
/// @param set DataSet
/// @param trait Index of trait to undo
/// @throws MethodException when negative value encountered
///
void LogTransform::UndoTrait(DataSet* set, int trait){
  int nInds = set->num_inds();
  double missingCoValue = set->get_missing_covalue();
  Sample * samp;
  for(int ind=0; ind < nInds; ind++){
    samp = set->get_sample(ind);
    if(samp->getTrait(trait)==missingCoValue)
      continue;
    samp->setTrait(Undo(samp->getTrait(trait)), trait);
  }
}


///
/// Transforms value and returns transformed value
/// @param val
/// @return transformed value
/// @throws MethodException when value is negative
///
double LogTransform::Transform(double val){
  if(val < 0)
    throw MethodException("Negative value passed to LogTransform Transform");
  return log(val+1);
}

///
/// Undoes transformation and returns original value
/// @param val
/// @return original value
/// @throws MethodException when value is negative
///
double LogTransform::Undo(double val){
  return exp(val)-1;
}


}
