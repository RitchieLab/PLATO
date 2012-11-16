//FisherExact.h

#include <math.h>
#include <vector>

using namespace std;

/// Calculates a 2X2 fisher's exact test
/// adapted from VassarStats javascript
/// http://faculty.vassar.edu/lowry/fisher.html

class FisherExact{

  public:
  
    /// constructor
    FisherExact();
    
    // calculate p-value based on the values in 2X2 table
    // a   b
    // c   d
    double fisher_2_2(unsigned a, unsigned b, unsigned c, unsigned d);
    
    // alternative call using a vector with a as index 0, b as index 1, etc.
    double fisher_2_2(vector<unsigned> & cells);

    
  private:
    // fills factorials so that no factorial must be calculated more than once
    void fill_factorials(unsigned n);
  
    vector<double> logfactorials;
    
};
