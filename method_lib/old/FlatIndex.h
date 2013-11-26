// FlatIndex.h

#ifndef FLATINDEX_H
#define FLATINDEX_H

#include <vector>
#include<iostream>
using namespace std;
namespace Methods{
///
/// Contains functionality for converting a multi-dimensional
/// array into a linear array.<br>
/// Adapted from Nate Barney's C MDR library
///


/// Converts multi-dimensional array into linear index
class FlatIndex{

  public:
    FlatIndex();
  
    unsigned int flatten_indexes(vector<int> & indexes);

    vector<unsigned int> decode_index(int index,  int numloci);
    vector<double> decode_index_double(int index,  int numloci);

    void set_genos_per_locus(unsigned int numGenos)
      {num_genos_per_locus = numGenos;}
      
    unsigned int get_genos_per_locus(){return num_genos_per_locus;}
    
    int get_size_array(unsigned int numLoci);
    
    vector<int> & get_valid_indexes(int numLoci){return includedIndexes.at(numLoci);}
       
    void set_included_indexes(int startLoc, int endLoc,
      bool includeAllCells, int missingValue);
    
    void set_non_missing_indexes(int comboSize, int missingValue);
    
    vector<vector<int> > get_included_indexes(){return includedIndexes;}
    
  private:
    unsigned int num_genos_per_locus;
    vector<vector<int> > includedIndexes;
};
};

#endif // FLATINDEX_H
