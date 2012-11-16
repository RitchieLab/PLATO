//LikelihoodRatio.h

#ifndef __LIKELIHOODRATIO_H__
#define __LIKELIHOODRATIO_H__

#include "DataSet.h"
#include "MethodException.h"
#include "ContingencyTable.h"
#include "StepOptions.h"

namespace Methods{
class LikelihoodRatio{

  public:
    LikelihoodRatio();
    LikelihoodRatio(DataSet* ds);
    
    /// runs logistic regression on the dataset set within the method
    void calculate(int locus);
    
    /// runs logistic regression on contingency table
    void calculate(ContingencyTable* table);
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);
    
    /// returns uncertainty coefficient
    float getLikelihoodRatio(){return likelihood_ratio;}
    
  private:
  
    void initialize();
    
    vector<Marker*> * markers;
    DataSet* dataset;
    float likelihood_ratio;
    ContingencyTable::TotalType total_type;
    
    map<string, ContingencyTable::TotalType> TotalTypeMap;

};
};

#endif
