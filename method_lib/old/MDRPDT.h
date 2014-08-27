// MDRPDT.h

#ifndef __MDRPDT_H__
#define __MDRPDT_H__

#include "DataSet.h"
#include "MethodException.h"
#include "StepOptions.h"

#include "mdrpdt/tstatistic.h"
#include "mdrpdt/ptestdistribution.h"
#include "mdrpdt/genotyperepository.h"
#include "mdrpdt/MDRPDTRepository.h"
#include "mdrpdt/evaluationmethod.h"
#include "mdrpdt/pdtmodel.h"

namespace Methods{

class MDRPDT{

  public:

    MDRPDT();
    MDRPDT(DataSet* ds);

    ~MDRPDT();

    /// runs MDR-PDT on the locus index passed
    void calculate(unsigned int locus);

    /// runs MDR-PDT on the loci passed in the vector
    void  calculate(vector<unsigned int>& loci);

    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);

    /// prepares DataSet for analysis by MDRPDT
    void resetDataSet(DataSet* ds){
      dataset = ds;
      if(params_set){
        Load(ds);
        if(num_ptests > 0)
          RunPTests();
      }
    }

    /// returns T statistics for training
    vector<float> getTstatTraining(){return t_stats_training;}

    /// returns T statistics for testing
    vector<float> getTstatTesting(){return t_stats_testing;}

    /// returns p value for the model
    float getPvalue(float mor_value){return dist->Evaluate(1,mor_value);}

    /// returns Matched Odds Ratios for the folds
    vector<float> getMatchedOddsRatio(){return mor_values;}

    /// returns average Matched Odds Ratio
    float getAvgMatchedOddsRatio();

  private:

    /// Loads dataset
    void Load(Methods::DataSet* ds);

    void initialize();

    int index_conversion(int loc);

    void runPDT(MdrPDT::PdtModel& model);

    DataSet* dataset;
    vector<float> t_stats_training, t_stats_testing, mor_values;
    float avg_mor;
    int num_ptests, random_seed;
    int xvCount;
    bool ptests_run, params_set;

    MdrPDT::MDRPDTRepository repo;
    MdrPDT::GenotypeRepository data;

  	MdrPDT::Evaluation::EvaluationMethod *eval;			///<Evaluation object (T-statistic)

  	MdrPDT::Distribution::PTestDistribution* dist;

  	void RunPTests();

};

}

#endif

