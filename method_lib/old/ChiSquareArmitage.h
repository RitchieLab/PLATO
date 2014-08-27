#ifndef CHISQUAREARMITAGE_H
#define CHISQUAREARMITAGE_H

#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

namespace Methods{
class ChiSquareArmitage{

  public:
    static double armitage(vector<vector<int> > & chiTotals);
};
};

#endif 
