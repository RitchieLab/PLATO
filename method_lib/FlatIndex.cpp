//FlatIndex.cpp

#include "FlatIndex.h"
namespace Methods{

//unsigned int FlatIndex::num_genos_per_locus=0;
//vector<vector<int> > FlatIndex::includedIndexes;
FlatIndex::FlatIndex(){
  num_genos_per_locus=0;
}


// Use: Convert multi-dimensional array indexes
//      into linear array index
// Arg: indexes -- vector containing index values
// Ret: index into linear array
unsigned int FlatIndex::flatten_indexes(vector<int> & indexes){
  
  unsigned int linear_index = 0;
  for(unsigned int i = 0; i < indexes.size(); i++){
    linear_index += indexes[i];
    
    if(i+1 < indexes.size())
      linear_index *= num_genos_per_locus;
  }

  return linear_index;
}



// Use: Converts the linear index value into
//      the vector (genotype) corresponding to that index
// Arg: index -- linear index for the genotype
//      numloci -- number of loci in genotype
// Ret: vector where each value is one locus
vector<unsigned int> FlatIndex::decode_index(int index, int numloci){
  vector<unsigned int> genotype(numloci,0);
  unsigned long locusIndex;
  for(int i=numloci; i; --i){
    locusIndex = i-1;
    genotype[locusIndex] = index % num_genos_per_locus;
    index /= num_genos_per_locus;
  }
  return genotype;
}

// Use: Converts the linear index value into
//      the vector (genotype) corresponding to that index
// Arg: index -- linear index for the genotype
//      numloci -- number of loci in genotype
// Ret: vector where each value is one locus
vector<double> FlatIndex::decode_index_double(int index,  int numloci){
  vector<double> genotype(numloci,0);
  unsigned long locusIndex;
  for(int i=numloci; i; --i){
    locusIndex = i-1;
    genotype[locusIndex] = index % num_genos_per_locus;
    index /= num_genos_per_locus;
  }
  return genotype;
}


// Use:  Returns total size of array needed
// Arg:  none
// Ret:  size
int FlatIndex::get_size_array(unsigned int numLoci){
  int size=1;
  for(unsigned int i=1; i <= numLoci; i++)
    size *= num_genos_per_locus; 
  return size;
}


// Use: Sets the included indexes.  These are the indexes
//      that are totalled when determining number of
//      individuals classified correctly.  When the option
//      to ignore 
// Arg: 
void FlatIndex::set_included_indexes(int startCombSize, int endCombSize,
  bool includeAllCells, int missingValue){

  // add empty vectors to the 2-D vector
  includedIndexes.assign(endCombSize+1, vector<int>(0,0));

  // determines if need to exclude any of the indexes
  // from being counted
  if(includeAllCells){
    int sizeVector, currCell;
    for(int currCombSize=startCombSize; currCombSize <= endCombSize;
      currCombSize++){
      sizeVector = FlatIndex::get_size_array(currCombSize);
      for(currCell=0; currCell<sizeVector; currCell++)
        includedIndexes[currCombSize].push_back(currCell);
    }
  }
  else{
    for(int currCombSize=startCombSize; currCombSize <= endCombSize;
      currCombSize++)
      set_non_missing_indexes(currCombSize, missingValue);
  } 
}


void FlatIndex::set_non_missing_indexes(int comboSize, int missingValue){
// cout << "start set_non_missing_indexes" << endl;
  int * lower_index = new int[comboSize];
  int * upper_index = new int[comboSize];

  for(int currLoc=0; currLoc<comboSize; currLoc++){
    lower_index[currLoc] = 0;
    upper_index[currLoc] = num_genos_per_locus;
  }
  
  int * indexes = new int[comboSize];
  int cur_depth = 0;
  int max_depth = comboSize;

  indexes[cur_depth] = lower_index[cur_depth];
  vector<int> loc_indexes; 
  bool containsMissing = false; 
  int i;
  
  while(1){
    if( indexes[cur_depth] < upper_index[cur_depth] ) {
      if( cur_depth == max_depth - 1 ) {
        // need to determine whether to include
        // or not here
        loc_indexes.clear();
        containsMissing = false;
        for(i=0; i<=cur_depth; i++){
          if(missingValue == indexes[i]){
            containsMissing = true;
            break;
          }
          loc_indexes.push_back(indexes[i]);
        }
        
        if(!containsMissing){
          includedIndexes[comboSize].push_back(FlatIndex::flatten_indexes(loc_indexes));
        }
        indexes[cur_depth]++;
      }
      else{
        cur_depth++;
        indexes[cur_depth]=lower_index[cur_depth];
      }
    }
    else{
      if(cur_depth > 0){
        indexes[--cur_depth]++;
      }
      else{
        break;
      }
    }
  }
  delete [] indexes;
  
  delete [] lower_index;
  delete [] upper_index;  
}

}
