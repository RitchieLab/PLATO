//NMI.cpp

#include "NMI.h"
// #include <iomanip>
namespace Methods{
///
/// Constructor
///
NMI::NMI(){
  dataset = NULL;
  initialize();
}


///
/// Alternative constructor
///
NMI::NMI(DataSet* ds){
  resetDataSet(ds);
  initialize();
}


///
/// Initialize starting variables
///
void NMI::initialize(){

  nmi_score = 0.0;
  transpose_on = false;
  
  total_type = ContingencyTable::Allele;
    
  TotalTypeMap["ADDITIVE"] = ContingencyTable::Additive;
  TotalTypeMap["RECESSIVE"] = ContingencyTable::Recessive;
  TotalTypeMap["DOMINANT"] = ContingencyTable::Dominant;
  TotalTypeMap["ALLELE"] = ContingencyTable::Allele;
  TotalTypeMap["GENOTYPE"] = ContingencyTable::Genotype;
  
  TransposeMap["TRUE"] = true; 
  TransposeMap["ON"] = true;
  TransposeMap["FALSE"] = false;
  TransposeMap["OFF"] = false;
}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void NMI::resetDataSet(DataSet* ds){
  dataset = ds;
  markers = dataset->get_markers();
}


///
/// Calculates odds ratio for the locus passed
/// @param locus 
///
void NMI::calculate(int locus){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  locus = (*markers).at(locus)->getLoc();
  
  ContingencyTable table;
  
  table.get_counts(locus, dataset);
  
  // this is opposite of expected as Will expects cases/ controls to be column
  // headers and I make them be the row identifiers
  // so when Will wants it to be NMIT then it doesn't need to be transposed
  if(!transpose_on){
    table = table.transpose();
  }
  calculate_score(table);
  
}


///
/// Calculates odds ratio
/// @param table Contingency table
///
void NMI::calculate(ContingencyTable* orig_table){

  ContingencyTable table; 
  // this is opposite of expected as Will expects cases/ controls to be column
  // headers and I make them be the row identifiers
  // so when Will wants it to be NMIT then it doesn't need to be transposed
  if(!transpose_on){
    ContingencyTable table = orig_table->transpose();
    calculate_score(table);
  }
  else
    calculate_score(*orig_table); 
}


///
/// Calculates and sets NMI score
/// @param table Contingency table to analyze
///
void NMI::calculate_score(ContingencyTable& table){

  table.set_current_totals(total_type);

  // for each row generate the entropy
  double H_B = 0.0;
  
  double H_B_A = 0.0;
  
  unsigned int num_rows = table.num_rows();
  unsigned int num_cols = table.num_cols();
  
  double total_in_table = table.get_total_in_table();

  unsigned int row, col;

  vector<vector<double> > table_float(num_rows, vector<double>(num_cols, 0.0));
  
  vector<double> row_totals(num_rows,0.0);
  vector<double> col_totals(num_cols,0.0);
  
  for(col=0; col<num_cols; col++){
    for(row=0; row<num_rows; row++){
      table_float.at(row).at(col) = table[row][col]/total_in_table;
      row_totals.at(row)+= table_float.at(row).at(col);
      col_totals.at(col)+= table_float.at(row).at(col);
    }
    H_B += (col_totals.at(col) * log2(col_totals.at(col)));
  }
  
  double row_value;
  
  for(row=0; row<num_rows; row++){
    row_value = 0.0;
    for(col=0; col<num_cols; col++){
      row_value += (double(table_float.at(row).at(col))/row_totals.at(row)) * log2(double(table_float.at(row).at(col))/row_totals.at(row));
    }

    H_B_A = H_B_A + row_totals.at(row) * row_value;
  }
  
  H_B_A = -H_B_A;
  H_B = -H_B;

  nmi_score = (H_B - H_B_A)/H_B;
}


///
/// Sets parameters for running method
/// @param options StepOptions
/// @throws MethodException if there is an invalid parameter
///
void NMI::set_parameters(StepOptions* options){

  if(TotalTypeMap.find(options->getNMITotalType()) != TotalTypeMap.end()){
    total_type = TotalTypeMap[options->getNMITotalType()];
  }
  else{
    throw MethodException(options->getNMITotalType() + 
      "= is not a recognized choice for NMI method");
  }
  transpose_on = options->getNMITransposed();
  
}
}
