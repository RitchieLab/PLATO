//UncertaintyCoefficient.cpp

#include "UncertaintyCoefficient.h"
namespace Methods{
///
/// Constructor
///
UncertaintyCoefficient::UncertaintyCoefficient(){
  dataset = NULL;
  initialize();
}


///
/// Alternative constructor
///
UncertaintyCoefficient::UncertaintyCoefficient(DataSet* ds){
  resetDataSet(ds);
  initialize();
}


///
/// Initialize starting variables
///
void UncertaintyCoefficient::initialize(){

//   maxLocusValue = 2;
  uncertainty_coeff = 0.0;
  
  total_type = ContingencyTable::Allele;
    
  TotalTypeMap["ADDITIVE"] = ContingencyTable::Additive;
  TotalTypeMap["RECESSIVE"] = ContingencyTable::Recessive;
  TotalTypeMap["DOMINANT"] = ContingencyTable::Dominant;
  TotalTypeMap["ALLELE"] = ContingencyTable::Allele;
  TotalTypeMap["GENOTYPE"] = ContingencyTable::Genotype;
  
//   missingValue = maxLocusValue+1;
//   LociComboLimit = 10;
//   LociComboMin = 1;
//   
//   set_threshold = 1.0;
  
}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void UncertaintyCoefficient::resetDataSet(DataSet* ds){
  dataset = ds;
//   calculate_set_threshold();
//   missingValue = dataset->get_missing_value();
  markers = dataset->get_markers();
//   setIndexConverter();
}


///
/// Calculates uncertainty coefficient for the locus passed
/// @param locus 
///
void UncertaintyCoefficient::calculate(int locus){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  locus = (*markers)[locus]->getLoc();
  
  ContingencyTable table;
  
  table.get_counts(locus, dataset);
  uncertainty_coeff = calculate_uc(&table);
}


void UncertaintyCoefficient::calculate(ContingencyTable* orig_table){
  uncertainty_coeff = calculate_uc(orig_table);
}


///
/// Calculates Uncertainty coefficient  
/// @param table Contingency table
///
float UncertaintyCoefficient::calculate_uc(ContingencyTable* orig_table){

  // use transposed table to get the Uncertainty Coefficient
  ContingencyTable table = orig_table->transpose();
  
  table.set_current_totals(total_type);

// table.output_current_grid();

  // calculate H(A) and H(B) (entropies of status and alleles)
  // then calculate T which is H(B)-H(B|a)
  // 2T(A,B) / (H(A) + H(B))
  
  unsigned int num_rows = table.num_rows();

  unsigned int num_cols = table.num_cols();
  double total_in_table = table.get_total_in_table();  
  unsigned int row, col;

  vector<vector<double> > table_float(num_rows, vector<double>(num_cols, 0.0));
  
  vector<double> row_totals(num_rows,0.0);
  vector<double> col_totals(num_cols,0.0);
  
  double H_B = 0.0, H_A = 0.0;
  
  // calculate row and column totals  
  for(row=0; row<num_rows; row++){
    for(col=0; col<num_cols; col++){
      table_float[row][col] = table[row][col]/total_in_table;
      row_totals[row]+= table_float[row][col];
      col_totals[col]+= table_float[row][col];
    }
    H_A += (row_totals[row] * log2(row_totals[row]));
  } 
  
  for(col=0; col<num_cols; col++){
    H_B += (col_totals[col] * log2(col_totals[col]));
  }
  
  H_B = -H_B;
  H_A = -H_A;

  double H_B_A = 0.0;
  // calculate H(B|a)
  
  double row_value;
  
  for(row=0; row<num_rows; row++){
    row_value = 0.0;
    for(col=0; col<num_cols; col++){
      row_value += (double(table_float[row][col])/row_totals[row]) * log2(double(table_float[row][col])/row_totals[row]);
    }

    H_B_A = H_B_A + row_totals[row] * row_value;
  }
  
  H_B_A = -H_B_A;
  
  double uc = 2 * (H_B - H_B_A) / (H_B + H_A);
  if(uc < 0){
    uc = 0;
  }
  return uc;
  
//   return 2 * (H_B - H_B_A) / (H_B + H_A);  
//   uncertainty_coeff = 2 * (H_B - H_B_A) / (H_B + H_A);
}

///
/// Sets parameters for running method
/// @param options StepOptions
/// @throws MethodException if there is an invalid parameter
///
void UncertaintyCoefficient::set_parameters(StepOptions* options){
  
  if(TotalTypeMap.find(options->getUncertaintyCoeffTotalType()) != TotalTypeMap.end()){
    total_type = TotalTypeMap[options->getUncertaintyCoeffTotalType()];
  }
  else{
    throw MethodException(options->getUncertaintyCoeffTotalType() + 
      " is not a recognized choice for Uncertainty Coefficient method");
  }
  
}
}
