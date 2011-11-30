#ifndef FAMILY_H
#define FAMILY_H

#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <cctype>
#include <algorithm>
#include "Options.h"
#include "Globals.h"
#include "Sample.h"
namespace Methods{

class Sample;

class Family{
	private:
		string famid;
		int famid_digit;
		string center;
		bool enabled;
		bool excluded;
		bool flag;
		bool flag_AF_CM;
		bool flag_AF_CF;
		bool has_parents;
		bool has_children;
		bool alphanumeric;
		bool me_error;
		int loc;
		int sysid;
		vector<Sample*> samples;
		vector<Sample*> founders;
		vector<Sample*> nonfounders;


		//For QFAM model
		double B;
		
	public:
		Family(){
			samples.resize(0);
			enabled = true;
			flag = false;
			has_parents = false;
			has_children = false;
			flag_AF_CM = false;
			flag_AF_CF = false;
			famid = "";
			center = "";
			loc = 0;
			famid_digit = -1;
			alphanumeric = true;
			excluded = false;
			me_error = false;
			B = 0;
		};

		~Family(){
		};
		
//		fam_data getFamStruct();

		string toString();

		//For QFAM model
		double getB(){return B;}
		void setB(double v){B = v;}
	
		void setSysid(int i){sysid = i;}
		int getSysid(){return sysid;}
		void setMeError(bool b){me_error = b;};
		bool hasMeError(){return me_error;};
		void setExcluded(bool b){excluded = b;};
		bool isExcluded(){return excluded;};
		bool isEnabled(){return enabled;};
		bool isFlagged(){return flag;};
		bool isFlaggedAFCM(){return flag_AF_CM;};
		bool isFlaggedAFCF(){return flag_AF_CF;};
		void setEnabled(bool v);
		void setFlag(bool v){flag = v;};
		void setFlagAFCM(bool v){flag_AF_CM = v;};
		void setFlagAFCF(bool v){flag_AF_CF = v;};
		string getCenter(){
			return center;
		};
		void setCenter(string cen){center = cen;};
		void addFounder(Sample* s);
		void addNonFounder(Sample* s);
		vector<Sample*>* getNonFounders(){return &nonfounders;};
		vector<Sample*>* getFounders(){return &founders;};
		int getLoc(){return loc;};
		void setLoc(int l){loc = l;};

		string getFamID();
		void AddInd(Sample*);

		void setFamID(string v){
			famid = v;
		};
		bool isAlphanumeric(){return alphanumeric;};
		void setAlphanumeric(bool b){alphanumeric = b;};
		string getFamIDOrig(){return famid;}
		void setFamDigit(int i){famid_digit = i;};
		int getFamID_digit(){
			if(famid_digit == -1){
				famid_digit = atoi(famid.c_str());
			}
			return famid_digit;
		};
		int getTotalInds(){return samples.size();};
		int getTotalEnabledInds();
		vector<Sample*>* getSamples(){return &samples;};
//		void resetGenderErrors();
};
};

#endif
