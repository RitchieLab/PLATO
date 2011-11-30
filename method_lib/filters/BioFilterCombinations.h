//BioFilterCombinations.h

#ifndef __BIOFILTERCOMBOS_H__
#define __BIOFILTERCOMBOS_H__

#include "DB.h"
#include "ResultSet.h"
#include <DataSet.h>

namespace Filters{

/// Get combinations from specified tables in database
class BioFilterCombinations{

  public:
    /// Contructor
    BioFilterCombinations();
    
    /// Alternative constructor
    BioFilterCombinations(string hostname, string username, string pass,
      string schema);
    
    /// Destructor
    ~BioFilterCombinations();
       
    /// Records information about set
    void fill_combinations(string tablename, uint combination_size);
    
    /// Retrieve combinations -- return false when all combinations returned
    bool get_combinations(vector<vector<string> >& combos);
        
  private:
  
    void initialize_combos(uint size);
    
    uint number_combos(string tablename);

    void add_combos();

    DB conn;
  
    vector<vector<string> > combinations;
    
    unsigned int combo_interval, combo_size, curr_res;
    
    string startsql, increment;

};

}

#endif
