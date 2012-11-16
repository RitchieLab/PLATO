//ResultSet.h

//Stores results and allows for addition
//and deletion of new results

#ifndef __RESULTSET_H__
#define __RESULTSET_H__

#include <list>
#include "DataDefinitions.h"

using namespace std;

  ///
  /// All results are stored as a combination of loci and   
  /// scores from each filter in the analysis
  ///
namespace Filters{
  /// Struct containing result information
  class Result{
    public:
    /// Combination of loci -- can be single locus
    LocusVec genoCombination;
    /// Scores at each index are for the corresponding analysis -- subfilter results are not stored
    vector<float> analysisScores;
    
    /// overload < for sort that sorts based on combination size
    int operator<(const Result& res) const{
    if(this->genoCombination.size() < res.genoCombination.size())
      return 1;
    else if(this->genoCombination.size() > res.genoCombination.size())
      return -1;
    else
      return 0;
    }
  };

  typedef list<Result> ResultSet;
  typedef list<Result>::iterator ResultIter;
}
#endif
