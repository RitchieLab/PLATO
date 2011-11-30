//SquareRootTransform.h

#ifndef __SQUAREROOTTRANSFORM_H__
#define __SQUAREROOTTRANSFORM_H__

#include "DataSet.h"
#include "StepOptions.h"
#include "MethodException.h"

///
/// Transforms continuous data by taking the square root of each covariate/trait.  Values
/// can be converted back to the original score by squaring each covariate.  Values
/// that match the designated missing value are ignored and all other values have 0.5 
/// added to them to avoid problem of zero values in the data.
/// Based on pg 415-417 in Biometry by Robert R. Sokal/F. James Rohlf (Third Edition)
///

namespace Methods{

class SquareRootTransform{

  public:
    /// Transforms all data in set.
    void TransformData(DataSet* set);
  
    /// Transforms a covariate
    void TransformCovar(DataSet* set, int covar);

    /// Transforms a trait
    void TransformTrait(DataSet* set, int trait);
  
    /// Returns all data back to original values
    void UndoTransform(DataSet* set);
  
    /// Undoes transformation of a covariate
    void UndoCovariate(DataSet* set, int covar);
  
    /// Undoes transformation of a trait
    void UndoTrait(DataSet* set, int trait);
    
  private:
    double Transform(double val);
    double Undo(double val);
  
};

}


#endif
