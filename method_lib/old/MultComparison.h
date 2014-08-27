/*
 * MultComparison.h
 *
 *  Created on: Feb 3, 2009
 *      Author: gilesjt
 */

#ifndef MULTCOMPARISON_H_
#define MULTCOMPARISON_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "Helpers.h"
#include "Options.h"
#include "StepOptions.h"

using namespace std;

namespace Methods{
	class Pair
{
public:

  double p;
  int l;

  bool operator< (const Pair & p2) const
  {
    return ( p < p2.p );
  }

};

class MultComparison{
public:
	MultComparison(StepOptions o){
		options = o;
	}
	~MultComparison(){};

	vector<double> get_genomic_control(){return pv_GC;}
	vector<double> get_sidak_single_step(){return pv_sidakSS;}
	vector<double> get_sidak_step_down(){return pv_sidakSD;}
	vector<double> get_fdr_bh(){return pv_BH;}
	vector<double> get_fdr_by(){return pv_BY;}
	vector<double> get_bonferroni(){return pv_bonf;}
	vector<double> get_holm(){return pv_holm;}

	double get_genomic_control(int i){
		if(pv_GC.size() == 0)
			return -1;
		return pv_GC[mapping[i]];
	}
	double get_sidak_single_step(int i){
		if(pv_sidakSS.size() == 0)
			return -1;
		return pv_sidakSS[mapping[i]];
	}
	double get_sidak_step_down(int i){
		if(pv_sidakSD.size() == 0)
			return -1;
		return pv_sidakSD[mapping[i]];}
	double get_fdr_bh(int i){
		if(pv_BH.size() == 0)
			return -1;
		return pv_BH[mapping[i]];}
	double get_fdr_by(int i){
		if(pv_BY.size() == 0)
			return -1;
		return pv_BY[mapping[i]];}
	double get_bonferroni(int i){
		if(pv_bonf.size() == 0)
			return -1;
		return pv_bonf[mapping[i]];}
	double get_holm(int i){
		if(pv_holm.size() == 0)
			return -1;
		return pv_holm[mapping[i]];}

	void calculate(vector<double>&, vector<int>&);
protected:
	  vector<double> pv_GC;
	  vector<double> pv_sidakSS;
	  vector<double> pv_sidakSD;
	  vector<double> pv_holm;
	  vector<double> pv_BH;
	  vector<double> pv_BY;
	  vector<double> pv_bonf;

	  map<int, int> mapping;

	  StepOptions options;
};
};
#endif /* MULTCOMPARISON_H_ */
