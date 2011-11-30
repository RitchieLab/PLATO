//LikelihoodRatio.cpp

#include "LikelihoodRatio.h"
namespace Methods{
///
/// Constructor
///
LikelihoodRatio::LikelihoodRatio(){
  dataset = NULL;
  initialize();
}


///
/// Alternative constructor
///
LikelihoodRatio::LikelihoodRatio(DataSet* ds){
  resetDataSet(ds);
  initialize();
}


///
/// Initialize starting variables
///
void LikelihoodRatio::initialize(){
  likelihood_ratio = 0.0;
  
  total_type = ContingencyTable::Allele;
    
  TotalTypeMap["ADDITIVE"] = ContingencyTable::Additive;
  TotalTypeMap["RECESSIVE"] = ContingencyTable::Recessive;
  TotalTypeMap["DOMINANT"] = ContingencyTable::Dominant;
  TotalTypeMap["ALLELE"] = ContingencyTable::Allele;
  TotalTypeMap["GENOTYPE"] = ContingencyTable::Genotype;
}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void LikelihoodRatio::resetDataSet(DataSet* ds){
  dataset = ds;
  markers = dataset->get_markers();
}


///
/// Calculates uncertainty coefficient for the locus passed
/// @param locus 
///
void LikelihoodRatio::calculate(int locus){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  locus = (*markers).at(locus)->getLoc();
  ContingencyTable table;
  
  table.get_counts(locus, dataset);
  calculate(&table);
  
}


///
/// Calculates likelihood ratio
/// @param table Contingency table
///
void LikelihoodRatio::calculate(ContingencyTable* table){
 
  table->set_current_totals(total_type);

  // only use with 2X2 table or 2X3 table
  double affected=0.0, unaffected=0.0;
  
  vector<double> totalInCell;
  for(unsigned int cell=0; cell<(*table)[0].size(); cell++){
    unaffected+=(*table)[0][cell];
    affected+=(*table)[1][cell];
    totalInCell.push_back((*table)[0][cell] + (*table)[1][cell]);
  }

  double totalinds = double(unaffected) + affected;
  double llr=0.0;
  // calculate expected for each cell
  for(unsigned int cell=0; cell<(*table)[0].size(); cell++){
    if( totalInCell.at(cell) > 0){
      llr += (*table)[0][cell] * log((*table)[0][cell] / (totalInCell.at(cell)*unaffected/totalinds));
      llr += (*table)[1][cell] * log((*table)[1][cell] / (totalInCell.at(cell)*affected/totalinds));
    }
  }

  likelihood_ratio = 2 * llr;
  
  if(likelihood_ratio < 0){
    likelihood_ratio = 0;
  }
  
}

///
/// Sets parameters for running method
/// @param options StepOptions
/// @throws MethodException if there is an invalid parameter
///
void LikelihoodRatio::set_parameters(StepOptions* options){
  
  if(TotalTypeMap.find(options->getLikelihoodRatioTotalType()) != TotalTypeMap.end()){
    total_type = TotalTypeMap[options->getLikelihoodRatioTotalType()];
  }
  else{
    throw MethodException(options->getLikelihoodRatioTotalType() + 
      "= is not a recognized choice for Uncertainty Coefficient method");
  }
  
}

}
