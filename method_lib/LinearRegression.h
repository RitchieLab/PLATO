//Liberally pulled from Plink and adapted to work in the Plato environment

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __LINEAR_H__
#define __LINEAR_H__

#include<vector>
//#include "plink.h"
#include "Marker.h"
#include "DataSet.h"
#include "MethodException.h"
#include "Model.h"

using namespace std;

namespace Methods{
class LinearRegression : public Model {

 public:

  LinearRegression();//Plink *);
  ~LinearRegression() { };

  void resetDataSet(DataSet* ds){data_set = ds;}

  void setDependent();

  void fitLM();
  void fitUnivariateLM();
  void noCluster();
  void setCluster(vector<int> & cl);
  void pruneY();

  void reset();
  vector<double> getCoefs();
  vector<double> getVar();
  vector<double> getPVals();
  vector<double> getPvalues(){return pvalues;}
  double getPValue();
  double getZ();
  double getTestCoef();
  int getCalcMissing();

  void calculate(int);
  void displayResults(ofstream &, Marker *);

  double calculateRSS();
  double calculateRSquared();
  double calculateAdjustedRSquared();
  double calculateMallowC(LinearRegression *);
  double calculateFTest(LinearRegression *);

  vector<double> getZs(){return ZS;}
  vector<string> getLabels(){return label;}

 private:

  vector<double> Y;
  vector<int> clst;
  vector<int> C;

  vector<double> se;
  vector<double> ZS;
  vector<double> pvalues;
  double chisq;

  vector<double> sig;
  vector<double> w;
  vector<vector<double> > u;
  vector<vector<double> > v;

  double varY;
  int nc;

  double RSS;

  //DataSet* data_set;

  bool cluster;
  void function(const int i, vector<double> & p );
  void setVariance();
  void display(vector<vector<double> > &);
  void display(vector<double>&);
  void display(vector<int>&);
};
};

#endif
