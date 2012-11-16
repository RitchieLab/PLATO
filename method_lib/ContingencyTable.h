 // ContingencyTable.h

#ifndef __CONTINGENCYTABLE_H__
#define __CONTINGENCYTABLE_H__

#include "MethodException.h"
#include "DataSet.h"

namespace Methods{
struct table_totals{
  vector<vector<float> > totals;
  float total_count;
};

///
/// Calculates contingency table based on
/// dataset.  Calculates both allele totals and genotype
/// totals and stores them.
///

/// Runs permutations for assigning significance to results
class ContingencyTable{

  public:
    /// Constructor
    ContingencyTable();
    
    /// Copy Constructor
    ContingencyTable(const ContingencyTable& origTable);

    /// Computes totals for the indicated locus from the dataset
    void get_counts(unsigned int curr_loc, DataSet* data);
       
    /// Returns value as set by the pointer
    vector<float>& operator[](unsigned int index){return current_totals->totals[index];}
    
    /// Sets type of results to return
    enum TotalType{
      Genotype,
      Allele,
      Dominant,
      Recessive,
      Additive
    };
    
    /// Set current type of totals to return using overloaded operator
    void set_current_totals(TotalType type);
    
    /// Returns number of rows in table
    unsigned int num_rows(){return current_totals->totals.size();}
    
    /// Returns number of columns  in table
    unsigned int num_cols(){return current_totals->totals[0].size();}
    
    /// Returns vector for use
    vector<vector<float> >& get_vector(){return current_totals->totals;}
    
    /// Transposes table and returns transposed one
    ContingencyTable transpose();
    
    /// Operator= overloaded
    ContingencyTable& operator=(const ContingencyTable& origTable);
    
    /// Returns total in table
    float get_total_in_table(){return current_totals->total_count;}
    
    /// Outputs current grid for checking values in table
    void output_current_grid(){output_grid(*current_totals, "grid");}
    
  private:
  
    void initialize();
    
    void copy(const ContingencyTable& origTable);
  
    void transpose_vector(vector<vector<float> >& orig,
      vector<vector<float> >& transposed);
  
    void output_grid(table_totals tot, string name);
  
    void correction(table_totals& tot);
  
    TotalType current_type;
  
    table_totals * current_totals;
    table_totals allele_totals;
    table_totals genotype_totals;
    table_totals dominant_totals;
    table_totals recessive_totals;
    table_totals additive_totals;
    unsigned int total_inds;
};
};

#endif
