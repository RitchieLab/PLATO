/* Taken and adapted from Plink*/

#ifndef __LINEAR_H__
#define __LINEAR_H__

#include<vector>
#include "plink.h"
#include "model.h"

using namespace std;

class LinearModel : public Model {

 public:

  LinearModel(Plink *);
  ~LinearModel() { };

  void setDependent();

  void fitLM();
  void fitUnivariateLM();
  void noCluster();
  void setCluster(vector<int> & cl);
  void pruneY();

  void reset();
  vector_t getCoefs();
  vector_t getVar();
  vector_t getPVals();
  double getPValue();

  void displayResults(ofstream &, Locus *);

  double calculateRSS(); 
  double calculateRSquared();
  double calculateAdjustedRSquared();
  double calculateMallowC(LinearModel *);
  double calculateFTest(LinearModel *);

 private:

  vector_t Y;
  vector<int> clst;
  vector<int> C;

  vector<double> se;
  double chisq;

  vector<double> sig;
  vector<double> w;
  vector<vector<double> > u;
  vector<vector<double> > v;

  double varY;
  int nc;

  double RSS;

  bool cluster;
  void function(const int i, vector<double> & p );  
  void setVariance();
};


#endif
