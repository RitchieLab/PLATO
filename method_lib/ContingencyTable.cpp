// ContingencyTable.cpp
#include "ContingencyTable.h"
namespace Methods{
///
/// Constructor
///
ContingencyTable::ContingencyTable(){
  initialize();
}

///
/// Initializes object
/// @return
///
void ContingencyTable::initialize(){
  set_current_totals(Genotype);
  total_inds = 0;
}


///
/// Copy Constructor
/// @param origTable Contingency Table to copy
///
ContingencyTable::ContingencyTable(const ContingencyTable& origTable){
  copy(origTable);
}


///
/// overloaded operator=
/// @param origTable ContingencyTable to copy
/// @return this object
///
ContingencyTable& ContingencyTable::operator=(const ContingencyTable& origTable){
  copy(origTable);
  return *this;
}


///
/// Performs copy
/// @param origTable Contingency Table to copy
/// @return
///
void ContingencyTable::copy(const ContingencyTable& origTable){
    allele_totals = origTable.allele_totals;
    genotype_totals = origTable.genotype_totals;
    dominant_totals = origTable.dominant_totals;
    recessive_totals = origTable.recessive_totals;
    set_current_totals(origTable.current_type);
}


///
/// Totals alleles and genotypes for use in filters
/// @param curr_loc index of locus being totaled
/// @param data DataSet to check
/// @return
///
void ContingencyTable::get_counts(unsigned int curr_loc, DataSet* data){
  
  unsigned int num_genotypes = data->get_max_locus()+1;
  unsigned int num_cells = num_genotypes;
  unsigned int num_alleles = data->get_max_allele();
    
  if(data->missing_data_present())
    ++num_cells;

  vector<vector<unsigned int> > counts(2, vector<unsigned int>(num_cells,0));
  
  genotype_totals.totals.assign(2, vector<float>(num_cells,0));
  allele_totals.totals.assign(2, vector<float>(num_alleles+1, 0));  
  
  dominant_totals = allele_totals;
  recessive_totals = allele_totals;

  int num_inds = data->num_inds();
  
  // count genotypes
  for(int curr_ind = 0; curr_ind < num_inds; ++curr_ind){
    ++counts[data->get_sample(curr_ind)->getAffected()][data->get_sample(curr_ind)->get_genotype(curr_loc)];
  }
  
  // remove last column if missing data is present so it will not be used
  if(data->missing_data_present()){
    genotype_totals.totals[0].pop_back();
    genotype_totals.totals[1].pop_back();
  }
  
  genotype_totals.total_count=0;
  
  // calculate total and set the totals in floats
  for(unsigned int i=0; i<=1; ++i){
    for(unsigned int j=0; j<=2; ++j){
      genotype_totals.totals[i][j]=counts[i][j];
      genotype_totals.total_count += genotype_totals.totals[i][j];
    }
  }

  // calculate allele totals
  // assume that have 3 genotypes per locus
  allele_totals.totals[0][0] = 2 * genotype_totals.totals[0][0] + genotype_totals.totals[0][1];
  allele_totals.totals[0][1] = 2 * genotype_totals.totals[0][2] + genotype_totals.totals[0][1];
  
  allele_totals.totals[1][0] = 2 * genotype_totals.totals[1][0] + genotype_totals.totals[1][1];
  allele_totals.totals[1][1] = 2 * genotype_totals.totals[1][2] + genotype_totals.totals[1][1];
  
  allele_totals.total_count = 2 * genotype_totals.total_count;
  
  int ref_allele_index = (*data->get_markers())[curr_loc]->getReferentIndex();
  dominant_totals.totals[0][0] = genotype_totals.totals[0][2*ref_allele_index]  + genotype_totals.totals[0][1];
  dominant_totals.totals[0][1] = genotype_totals.totals[0][2-2*ref_allele_index];
  dominant_totals.totals[1][0] = genotype_totals.totals[1][2*ref_allele_index] + genotype_totals.totals[1][1];
  dominant_totals.totals[1][1] = genotype_totals.totals[1][2-2*ref_allele_index];
  
//   //calculate dominant model totals
  dominant_totals.total_count = genotype_totals.total_count;
 
  // calculate recessive model totals
  recessive_totals.totals[0][0] = genotype_totals.totals[0][2*ref_allele_index]; 
  recessive_totals.totals[0][1] = genotype_totals.totals[0][2-2*ref_allele_index] + genotype_totals.totals[0][1];
  
  recessive_totals.totals[1][0] = genotype_totals.totals[1][2*ref_allele_index];
  recessive_totals.totals[1][1] = genotype_totals.totals[1][2-2*ref_allele_index] + genotype_totals.totals[1][1];  
  recessive_totals.total_count = genotype_totals.total_count;


  current_totals = &allele_totals;


  // Do correction for all totals by adding 0.5 to each cell
  correction(allele_totals);
  correction(genotype_totals);
  correction(dominant_totals);
  correction(recessive_totals);
  
  additive_totals = genotype_totals;
  
  if(ref_allele_index == 0){
    additive_totals.totals[0][0] = genotype_totals.totals[0][2];
    additive_totals.totals[0][1] = genotype_totals.totals[0][1];
    additive_totals.totals[0][2] = genotype_totals.totals[0][0];
    
    additive_totals.totals[1][0] = genotype_totals.totals[1][2];
    additive_totals.totals[1][1] = genotype_totals.totals[1][1];
    additive_totals.totals[1][2] = genotype_totals.totals[1][0];
    additive_totals.total_count = genotype_totals.total_count;    
  }
  
}


///
/// Outputs grid
///
void ContingencyTable::output_grid(table_totals tot, string name){
  cout << name << endl;
  for(unsigned int i=0; i<tot.totals.size(); i++){
    for(unsigned int j=0; j<tot.totals[i].size(); j++){
      cout << tot.totals[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

///
/// Increments each cell by 0.5 for the correction
/// @param tot table_totals
///
void ContingencyTable::correction(table_totals& tot){
  unsigned int num_cols = tot.totals[0].size();
  unsigned int j;
  for(unsigned int i=0; i<tot.totals.size(); ++i){
    for(j=0; j<num_cols; ++j){
      tot.totals[i][j] += 0.5;
      tot.total_count += 0.5;
    }
  }
}


///
/// Transposes current table and returns new
/// @return Contingency table with results transposed
///
ContingencyTable ContingencyTable::transpose(){
  
  ContingencyTable new_table;
  
  // convert columns into rows
  transpose_vector(genotype_totals.totals, new_table.genotype_totals.totals);
  transpose_vector(recessive_totals.totals, new_table.recessive_totals.totals);
  transpose_vector(dominant_totals.totals, new_table.dominant_totals.totals);
  transpose_vector(allele_totals.totals, new_table.allele_totals.totals);

  new_table.allele_totals.total_count = allele_totals.total_count;
  new_table.recessive_totals.total_count = recessive_totals.total_count;
  new_table.genotype_totals.total_count = genotype_totals.total_count;
  new_table.dominant_totals.total_count = dominant_totals.total_count;
  new_table.additive_totals = new_table.genotype_totals;

  new_table.set_current_totals(current_type);

  return new_table;
}



///
/// Converts original vector into transposed one
/// @param orig original vector
/// @param transposed vector with transposed 2-D vector
///
void ContingencyTable::transpose_vector(vector<vector<float> >& orig,
  vector<vector<float> >& transposed){
  
  unsigned int num_rows = orig.size();
  unsigned int num_cols = orig[0].size();
  unsigned int orig_col;
  
  transposed.assign(num_cols, vector<float>(0,0));
  
  for(unsigned int orig_row=0; orig_row<num_rows; ++orig_row){
    for(orig_col=0; orig_col<num_cols; ++orig_col){
      transposed[orig_col].push_back(orig[orig_row][orig_col]);
    }
  }
  
}


///
/// Set current type of totals to return using overloaded operator
/// @param type TotalType
/// @return
/// @throws MethodExcept when 
///
void ContingencyTable::set_current_totals(TotalType type){
  current_type = type;
  
  switch(type){
    case Allele:
      current_totals = &allele_totals;
      break;
    case Genotype:
      current_totals = &genotype_totals;
      break;
    case Dominant:
      current_totals = &dominant_totals;
      break;
    case Recessive:
      current_totals = &recessive_totals;
      break;
    case Additive:
      current_totals = &genotype_totals;
      break;
    default:
      throw MethodException("Not valid TotalType");
  };
}


}
