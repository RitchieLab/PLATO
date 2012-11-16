//ChiSquare.cpp

#include "ChiSquareAllelic.h"

///
///probability of chi sqaure value 
///ALGORITHM Compute probability of chi square value.
///  Adapted from:
///    Hill, I. D. and Pike, M. C.  Algorithm 299
///    Collected Algorithms for the CACM 1967 p. 243
///  Updated for rounding errors based on remark in
///    ACM TOMS June 1985, page 185
///@param x obtained chi-square value
///@param df degrees of freedom
///@return p value
///
double ChiSquareAllelic::pfromchi(double x, int df){
  double  a, s, y=0;
  double  e, c, z;
  int   even;     /* true if df is an even number */
  
  if (x <= 0.0 || df < 1)
    return (1.0);
  
  a = 0.5 * x;
  even = (2*(df/2)) == df;
  if (df > 1)
    y = ex (-a);
  s = (even ? y : (2.0 * poz (-sqrt (x))));
  if (df > 2){
    x = 0.5 * (df - 1.0);
    z = (even ? 1.0 : 0.5);
    if (a > BIGX)
      {
      e = (even ? 0.0 : LOG_SQRT_PI);
      c = log (a);
      while (z <= x)
        {
        e = log (z) + e;
        s += ex (c*z-a-e);
        z += 1.0;
        }
      return (s);
      }
    else
      {
      e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
      c = 0.0;
      while (z <= x)
        {
        e = e * (a / z);
        c = c + e;
        z += 1.0;
        }
      return (c * y + s);
      }
    }
  else
    return (s);
}

///
/// probability of normal z value
/// Adapted from a polynomial approximation in:
///  Ibbetson D, Algorithm 209
///  Collected Algorithms of the CACM 1963 p. 616
///  Note:
///    This routine has six digit accuracy, so it is only useful for absolute
///    z values < 6.  For z values >= to 6.0, poz() returns 0.0.
/// @param  poz Normal z value
/// @return  probability
///
double ChiSquareAllelic::poz (double z)
{
  double x = 0.0;
  
  if (z != 0.0)
  {
    double y = 0.5 * fabs (z);
    if (y >= 3.0)
      x = 1.0;
    else if (y < 1.0) {
      double w = y*y;
      x = ((((((((0.000124818987 * w
        -0.001075204047) * w +0.005198775019) * w
        -0.019198292004) * w +0.059054035642) * w
        -0.151968751364) * w +0.319152932694) * w
        -0.531923007300) * w +0.797884560593) * y * 2.0;
    }
    else {
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y
        +0.000152529290) * y -0.000019538132) * y
        -0.000676904986) * y +0.001390604284) * y
        -0.000794620820) * y -0.002034254874) * y
        +0.006549791214) * y -0.010557625006) * y
        +0.011630447319) * y -0.009279453341) * y
        +0.005353579108) * y -0.002141268741) * y
        +0.000535310849) * y +0.999936657524;
    }
  }
  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
} // poz


///
/// Calculate chi square <br>
/// Assumes binary status with each individual being either a 0 or 1
/// @param chiTotals two dimensional vector with each row being a status and
///        each column a genotype.  Values represent totals for
///        that status and genotype
/// @param df pointer to integer that will hold degrees of freedom for this
///        model -- it is an output argument
/// Ret: chi square
///
vector<vector<double> > ChiSquareAllelic::expecteds(vector<vector<int> > & chiTotals){
  int affected=0, unaffected=0;
  vector<int> totalInCell;
  for(unsigned int cell=0; cell<chiTotals[0].size(); cell++){
    unaffected+=chiTotals[0][cell];
    affected+=chiTotals[1][cell];
    totalInCell.push_back(chiTotals[0][cell] + chiTotals[1][cell]);
  }

  float totalinds = float(unaffected) + affected;
  double chiscore=0.0;
  // calculate expected for each cell
  vector<vector<double> > expected;
  expected.resize(2);
  for(unsigned int cell=0; cell<chiTotals[0].size(); cell++){
    if( totalInCell[cell] > 0){
      float unaffexpected = float(totalInCell[cell]) * unaffected / totalinds;
	  expected[0].push_back(unaffexpected);
      float affexpected = float(totalInCell[cell]) * affected / totalinds;
	  expected[1].push_back(affexpected);
	}
  }
	
  return expected;
}

double ChiSquareAllelic::chisquare(vector<vector<int> > & chiTotals){
  int affected=0, unaffected=0;
  vector<int> totalInCell;
  for(unsigned int cell=0; cell<chiTotals[0].size(); cell++){
    unaffected+=chiTotals[0][cell];
    affected+=chiTotals[1][cell];
    totalInCell.push_back(chiTotals[0][cell] + chiTotals[1][cell]);
  }

  float totalinds = float(unaffected) + affected;
  double chiscore=0.0;
  // calculate expected for each cell
  for(unsigned int cell=0; cell<chiTotals[0].size(); cell++){
    if( totalInCell[cell] > 0){
      float unaffexpected = float(totalInCell[cell]) * unaffected / totalinds;
      float affexpected = float(totalInCell[cell]) * affected / totalinds;
      chiscore += pow((chiTotals[0][cell]-unaffexpected),2) / unaffexpected;
      chiscore += pow((chiTotals[1][cell]-affexpected),2) / affexpected;
    }
  }
  return chiscore;
}

///
/// Calculates degrees of freedom <br>
/// assumes binary status so only 2 rows and
/// only need to calculate number of columns in the table  
/// @param  numloci number of loci in model
/// @param  numLociValues number of possible loci values
/// @return degrees of freedom
///
int ChiSquareAllelic::calcdf(int numloci, int numLociValues){
  return int(pow(double(numLociValues), numloci)) - 1;
}
