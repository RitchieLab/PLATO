//OddsRatio.cpp

#include "OddsRatio.h"
namespace Methods{
///
/// Constructor
///
OddsRatio::OddsRatio(){
  dataset = NULL;
  initialize();
}

OddsRatio::~OddsRatio(){

}


///
/// Alternative constructor
///
OddsRatio::OddsRatio(DataSet* ds){
  resetDataSet(ds);
  initialize();
}


///
/// Initialize starting variables
///
void OddsRatio::initialize(){

  odds_ratio = 0.0;

  total_type = ContingencyTable::Allele;
  TotalTypeMap["ALLELE"] = ContingencyTable::Allele;
}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void OddsRatio::resetDataSet(DataSet* ds){
  dataset = ds;
  markers = dataset->get_markers();
}


///
/// Calculates odds ratio for the locus passed
/// @param locus
///
void OddsRatio::calculate(int locus){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  locus = (*markers).at(locus)->getLoc();

  ContingencyTable table;

  table.get_counts(locus, dataset);
  calculate(&table);

}


///
/// Calculates odds ratio
/// @param table Contingency table
///
void OddsRatio::calculate(ContingencyTable* table){

  table->set_current_totals(total_type);

  //AD/BC
  // only use with 2X2 table

  float numerator, denominator;

  if(((*table)[0][0] + (*table)[1][1]) >= ((*table)[0][1] + (*table)[1][0])){
    numerator = (*table)[0][0] * (*table)[1][1];
    denominator = (*table)[0][1] * (*table)[1][0];
  }
  else{
    numerator = (*table)[0][1] * (*table)[1][0];
    denominator = (*table)[0][0] * (*table)[1][1];
  }

  odds_ratio = numerator / denominator;

}

///
/// Sets parameters for running method
/// @param options StepOptions
/// @throws MethodException if there is an invalid parameter
///
void OddsRatio::set_parameters(StepOptions* options){

  if(TotalTypeMap.find(options->getOddsRatioTotalType()) != TotalTypeMap.end()){
    total_type = TotalTypeMap[options->getOddsRatioTotalType()];
  }
  else{
    throw MethodException(options->getOddsRatioTotalType() +
      "= is not a recognized choice for Uncertainty Coefficient method");
  }

}
}
