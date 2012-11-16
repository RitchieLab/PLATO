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


#ifndef __MODEL_H__
#define __MODEL_H__

#include<vector>
#include <set>
#include "StepOptions.h"
#include "DataSet.h"
#include "MethodException.h"
//#include "plink.h"

using namespace std;


namespace Methods{
class int2{
public:
	int p1;
	int p2;
	int2(){p1 = p2 = 0;}
	int2(int a, int b){p1 = a; p2= b;}

	bool operator< (const int2 &b) const{
		return (p1 < b.p1 || (p1 == b.p1 && p2 < b.p2) );
	}
	bool operator== (const int2 &b) const{
		return (p1 == b.p1 && p2 == b.p2);
	}
};


class Model {

 public:

  Model();
  virtual ~Model() { };

  virtual void setDependent() = 0;
  virtual void pruneY() = 0;
  virtual void fitLM() = 0;
  virtual vector<double> getCoefs() = 0;
  virtual vector<double> getVar() = 0;
  virtual vector<double> getPVals() = 0;
  virtual void displayResults(ofstream &, Marker *) = 0;
  virtual void fitUnivariateLM() = 0;


  void setMissing();
  vector<bool> getMissing();
  void setMissing(vector<bool>&);
  void yokeMissing(Model *);
  void setHaploid();
  void setX();
  void setDominant();
  void setRecessive();
  void hasSNPs(bool);
  void addAdditiveSNP(int m){addAdditiveSNP(data_set->get_locus(m));};
  void addAdditiveSNP(Marker*);
  void addDominanceSNP(int m){addDominanceSNP(data_set->get_locus(m));};
  void addDominanceSNP(Marker*);
  void addHaplotypeDosage(set<int>&);
  void addSexEffect();
  bool isSexInModel();
  void addCovariate(int);
  void addInteraction(int,int);
  void buildDesignMatrix();
  bool checkVIF();
  vector<bool> validParameters();
  bool isValid() { return all_valid; }
  double getStatistic();
  //  double getPValue();
  double linearHypothesis(vector<vector<double> > &, vector<double> &);
  int Ysize() { return nind; }
  int getNP() { return np; }
  void setValid() { all_valid = true; }
  vector<string> label;
  vector<int> order;
  vector<int> type;
  int testParameter;
  void setThreshold(string s){
	  options.setUp(s);
	  setOptions(options);
  }
  void setOptions(StepOptions o){
	  options = o;
	  if(options.getLinRConditionFile().size() > 0){
		  options.readLinRConditionFile(data_set->get_markers());
	  }
	  if(options.getLinRConditionString().size() > 0){
		  options.parseLinRConditionList(data_set->get_markers());
	  }
  }
  void resetDataSet(DataSet* ds){data_set = ds;}

  void nullCoefs();

  // Independent variables (can be directly manipulated...)
  vector<vector<double> > X;

  void setIsNull(bool b){isnull = b;}
  bool isNull(){return isnull;}

 protected:

  //Plink * P;
	 StepOptions options;
  // Missing flag
  vector<bool> miss;

  int nind;
  int np;  // Main effects + interaction + intercept

  bool has_snps;

  vector<bool> xchr;
  vector<bool> haploid;

  bool sex_effect;

  vector<bool> valid;
  bool all_valid;

  vector<double> coef;  // beta
  vector<vector<double> > S;     // Sigma


  // Term types

  enum terms { INTERCEPT,
	       ADDITIVE,
	       DOMDEV,
	       HAPLOTYPE,
	       SEX,
	       COVARIATE,
	       INTERACTION,
               QFAM };


  double buildIntercept();
  double buildAdditive(Sample *, int);
  double buildDominance(Sample *, int);
  double buildHaplotype(int, int);
  double buildSex(Sample *);
  double buildCovariate(Sample *, int);
  double buildInteraction(Sample *, int, vector<double> &);
  double buildQFAM(Sample *);

  bool skip;

  // List of additive SNP effects
  // assuming SNP major mode

  vector<Marker*> additive;

  int mAA;
  int mAB;
  int mBB;

  double mA, mB;

  // List of dominance deviation SNP effects

  vector<Marker*> dominance;

  // List of covariates (clist)

  vector<int> covariate;

  // List of pairwise interactions
  // ( indexing previously specified components, 1,2,..)

  vector<int2> interaction;

  // List of sets of haplotypes

  vector<set<int> > haplotype;

  DataSet* data_set;

  void display(vector<vector<double> > &);
  void display(vector<double> &);
  void display(vector<int> &);


  bool isnull;
};
};


#endif
