//UncertaintyCoefficient.h

#ifndef __UNCERTAINTYCOEFF_H__
#define __UNCERTAINTYCOEFF_H__

#include "DataSet.h"
#include "MethodException.h"
#include "ContingencyTable.h"
#include "StepOptions.h"

namespace Methods{
class UncertaintyCoefficient{

  public:
    UncertaintyCoefficient();
    UncertaintyCoefficient(DataSet* ds);
    
    /// runs logistic regression on the dataset set within the method
    void calculate(int locus);
    
    /// runs logistic regression on contingency table
    void calculate(ContingencyTable* table);
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);
    
    /// returns uncertainty coefficient
    float getUncertaintyCoefficient(){return uncertainty_coeff;}
    
  private:
    void initialize();
    
    float calculate_uc(ContingencyTable* orig_table);
    
    vector<Marker*> * markers;
    DataSet* dataset;
    float uncertainty_coeff;
    ContingencyTable::TotalType total_type;
    
    map<string, ContingencyTable::TotalType> TotalTypeMap;

};

};
#endif
