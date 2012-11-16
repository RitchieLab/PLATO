//FisherExact.cpp
#include "FisherExact.h"
namespace Methods{
///
/// Initializes logfactorials vector
///
FisherExact::FisherExact(){
  // set first 2 factorials to be zero
  logfactorials.push_back(0);
  logfactorials.push_back(log(double(1)));
}


///
/// Caches factorial results so that on repeated
/// calls to calculate Fisher's Exact test
/// these do not have to be calculated repeatedly
/// @param n maximum factorial to calculate
/// @return
///
void FisherExact::fill_factorials(unsigned n){

  // only add the factorials that have not
  // already been calculated
  for(unsigned i=logfactorials.size(); i<=n; i++){
    logfactorials.push_back(log(double(i)) + logfactorials[i-1]);
  }

}


///
/// Calculates and returns the two-tailed p-value
/// for this 2X2 contingency table as below
/// a  b
/// c  d
/// @param a
/// @param b
/// @param c
/// @param d
/// @return the two-tailed p-value
///
double FisherExact::fisher_2_2(unsigned a, unsigned b,
  unsigned c, unsigned d){

  double p_two = 0;
  int ab = a+b;
  int cd = c+d;
  int ac = a+c;
  int bd = b+d;
  int n = ab+cd;
  int np = n+1;

  // check to see if need to fill more factorials
  // in array
  fill_factorials(n);

  double disc = fabs((a/double(ab))-(c/double(cd)));
  disc = floor(disc*1000000);

  double factn;

  factn = logfactorials[n];
  double num = logfactorials[ab]+logfactorials[cd]+
    logfactorials[ac]+logfactorials[bd];

  int ax,bx,cx,dx;
  double denom;
  // go through for every data point
  // calculate probability when the degree of disproportion is
  // equal or greater than that of the observed array
  for (int j = 0; j < np; j++) {
    ax = j;
    if(ax>int(ab) || ax>int(ac))
      continue;
    bx = ab-ax;
    if (bx<0) continue;
    cx = ac-ax;
    if (cx<0) continue;
    dx = cd-cx;
    if (dx<0) continue;
    double newdisc = fabs((ax/double(ab))-(cx/double(cd)));
    newdisc = floor(newdisc*1000000);
    if (newdisc< disc) continue;

    denom = factn + logfactorials[ax] + logfactorials[bx] +
      logfactorials[cx] + logfactorials[dx];
    p_two += exp(num-denom);
  }
  
  if (p_two>.99999999) {p_two = 1.0;}

  return p_two;     
}

///
/// alternative call using a vector 
/// with a as index 0, b as index 1, etc.
/// @param cells vector of totals for table
/// @return two tailed p-valuef
///
double FisherExact::fisher_2_2(vector<unsigned> & cells){
  return fisher_2_2(cells[0], cells[1], cells[2], cells[3]);
}

}
