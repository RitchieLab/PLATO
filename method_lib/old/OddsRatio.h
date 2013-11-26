//OddsRatio.h

#ifndef __ODDSRATIO_H__
#define __ODDSRATIO_H__

#include "DataSet.h"
#include "MethodException.h"
#include "ContingencyTable.h"
#include "StepOptions.h"

namespace Methods{
class OddsRatio{

  public:
    OddsRatio();
    OddsRatio(DataSet* ds);
    virtual ~OddsRatio();

    /// runs logistic regression on the dataset set within the method
    void calculate(int locus);

    /// runs logistic regression on contingency table
    void calculate(ContingencyTable* table);

    /// sets DataSet
    void resetDataSet(DataSet* ds);

    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);

    /// returns uncertainty coefficient
    float getOddsRatio(){return odds_ratio;}

    /// returns 1 when first score is better than second one
    virtual int score_better(float first_score, float second_score){
      return first_score>second_score?1:0;}

  private:

    void initialize();

    vector<Marker*> * markers;
    DataSet* dataset;
    float odds_ratio;
    ContingencyTable::TotalType total_type;

    map<string, ContingencyTable::TotalType> TotalTypeMap;

};

};
#endif
