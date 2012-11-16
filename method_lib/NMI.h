//NMI.h

#ifndef __NMICALC_H__
#define __NMICALC_H__

#include "DataSet.h"
#include "MethodException.h"
#include "ContingencyTable.h"
#include "StepOptions.h"

namespace Methods{
class NMI{

  public:
    NMI();
    NMI(DataSet* ds);
    
    /// runs logistic regression on the dataset set within the method
    void calculate(int locus);
    
    /// runs logistic regression on contingency table
    void calculate(ContingencyTable* table);
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);
    
    /// returns uncertainty coefficient
    float getNMI(){return nmi_score;}
    
  private:
  
    void initialize();
    
    void calculate_score(ContingencyTable& table);
    
    vector<Marker*> * markers;
    DataSet* dataset;
    float nmi_score;
    ContingencyTable::TotalType total_type;
    
    bool transpose_on;
    map<string, ContingencyTable::TotalType> TotalTypeMap;
    map<string, bool> TransposeMap;
    
};

};
#endif
