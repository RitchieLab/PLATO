//MDR.h

#ifndef __MDR_H__
#define __MDR_H__

#include "DataSet.h"
#include "FlatIndex.h"
#include "StepOptions.h"
#include "MethodException.h"

///
/// Calculates MDR in method library
/// used by plato and wasp.  This method currently
/// only accepts snps for markers.
///

/// Method class that performs MDR
namespace Methods{
class MDR{

  public:
    MDR();
    MDR(DataSet* ds);
  
    /// run MDR on the indicated loci
    void calculate(vector<unsigned int>& loci, DataSet* set);
  
    /// runs MDR on the dataset set within the method
    void calculate(vector<unsigned int> loci);
  
    /// runs MDR on mixture of SNPs, covariates and traits
    void calculate(vector<unsigned int> loci, 
      vector<unsigned int>& covars, vector<unsigned int> & traits);  
  
    /// returns the balanced accuracy of the model
    double getBalancedAcc(){return mod.balanced_accuracy;}
    
    /// returns the number of correctly identified affected individuals
    int getTruePositives(){return mod.classhigh;}
    
    /// returns the number of correctly identified affected individuals
    int getFalsePositives(){return mod.misclasshigh;}
    
    /// returns the number of correctly identified affected individuals
    int getTrueNegatives(){return mod.classlow;}

    /// returns the number of correctly identified affected individuals
    int getFalseNegatives(){return mod.misclasslow;}

    /// returns the number of affected in each cell of model (as linear vector)
    vector<unsigned int> getAffectedCellTotals(){return mod.affected;}
    
    /// returns the number of unaffected in each cell of model (as linear vector)
    vector<unsigned int> getUnaffectedCellTotals(){return mod.unaffected;}
  
    /// set parameters for method using StepOptions class
    void set_parameters(StepOptions* options);
    
    /// sets whether to calculate the threshold for each model evaluated
    void setCalcThreshModel(bool calc){calc_each_threshold=calc;}
    
    /// sets DataSet 
    void resetDataSet(DataSet* ds);
    
    /// Returns the index into the affected and unaffected totals vectors for the genotypes passed
    unsigned int get_totals_index(vector<int>& genotypes){
      return indexConverter.flatten_indexes(genotypes);
    }
    
    /// Returns a vector containing the genotype values for a specific index in the totals vectors
    vector<unsigned int> genotypes_at_index(int index,  int numloci){
      return indexConverter.decode_index(index, numloci);
    }
    
  private:
    
    struct Model{
      vector<unsigned int> affected, unaffected;
      int misclasshigh, misclasslow, classhigh, classlow, totaltiecells;
      float threshold, balanced_accuracy;
    };
  
    DataSet* dataset;
  
    void initialize();
    void runMDR(unsigned int model_size);
    void calculate_set_threshold();
    float calculate_model_threshold(vector<unsigned int> & lociComb);
    float calculate_model_threshold(vector<unsigned int> loci, 
      vector<unsigned int>& covars, vector<unsigned int> & traits);

    void distributeInds(vector<unsigned int> & lociComb);
    
    void distributeInds(vector<unsigned int>& loci, 
      vector<unsigned int>& covars, vector<unsigned int> & traits, int miss);
    
    void calculateStats(unsigned int combSize);
    void setIndexConverter(int max_locus, bool any_missing_data, int missing_val);
    int set_converter_covars(vector<unsigned int> & covars,
      vector<unsigned int> & traits);
    void convert_loci_indexes(vector<unsigned int> & loci);
    
    int maxLocusValue, tieCellValue, missingCoValue;
    
    unsigned int missingValue, LociComboMin, LociComboLimit;
    bool calc_each_threshold;
    vector<vector< int> > includedIndexes;
    FlatIndex indexConverter;
    Model mod;
    vector<Marker*> * markers;
    
    float set_threshold;

};
};
#endif
