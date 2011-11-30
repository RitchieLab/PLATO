#ifndef CHISQUAREALLELIC_H
#define CHISQUAREALLELIC_H

#include <vector>
#include <math.h>
#include <iostream>

#define  CHI_EPSILON     0.000001    /* accuracy of critchi approximation */
#define  CHI_MAX     99999.0         /* maximum chi square value */

#define  LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define  I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define  BIGX           20.0         /* max value to represent exp (x) */
#define  ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))

using namespace std;

namespace Methods{
///
/// Class for calculating chi square and 
/// p values. It is used by ChiSquareAllelicFilter
///

/// Used by ChiSquareAllelicFilter
class ChiSquareAllelic{

  public:
    static double pfromchi(double x, int df);
    static vector<vector<double> >  expecteds(vector<vector<int> > & chiTotals);
    static double chisquare(vector<vector<int> > & chiTotals);
    static int calcdf(int numloci, int numLociValues);
    
  private:
    static double poz (double z);

};
};

#endif 
