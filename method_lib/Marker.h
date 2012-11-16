/*include for Marker.cc*/

#ifndef MARKER_H
#define MARKER_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include <occi.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <bitset>
#include "Globals.h"
#include "General.h"
#include "Options.h"
//#include "Helper.h"
//using namespace oracle::occi;
using namespace std;


namespace Methods{
class Marker{
	private:
		string enzyme;
		string probe_id;
		string allele1;
		string allele2;
		//list of alleles
		vector<string> alleles;
		float maf;
		bool enabled;
		bool flag;
		bool freqflag;
		int loc;
		int sysprobe;
		//micro satellite?
		bitset<1> micro_sat;
		double cM;
		//additional snp info stored here
		map<string,string> details;
		// referent allele for use in models -- represents disease allele
		string referent_allele;
		int ref_all_index;

	public:
		int chrom;
		int bploc;
		string rsid;
		Marker(string chr, string probe, int bp){
			//sysprobe = s;
			probe_id = probe;
			if(chr == "X" || chr == "x"){
				chrom = opts::_CHRX_;
			}
			else if(chr == "Y" || chr == "y"){
				chrom = opts::_CHRY_;
			}
			else if(chr == "XY" || chr == "xY" || chr == "Xy" || chr == "xy"){
				chrom = opts::_CHRXY_;
			}
			else{
				chrom = atoi(chr.c_str());
			}
			bploc = bp;
			enzyme = "";
			rsid = probe;
			allele1 = "";
			allele2 = "";
			enabled = true;
			flag = false;
			loc = -1;
			maf = -1;
			freqflag = false;
  		ref_all_index = -1;
			referent_allele = "";

		};
		Marker(){
			probe_id = "";
			rsid = "";
			chrom = -1;
			enabled = false;
			flag = false;
			enzyme = "";
			allele1 = "";
			allele2 = "";
			bploc = -1;
			loc = -1;
			freqflag = false;
			maf = -1;
			ref_all_index = -1;
			referent_allele = "";

		};
		~Marker(){};

		string toString();

		void setSysprobe(int i){sysprobe = i;};
		int getSysprobe(){return sysprobe;};
		void setLoc(int i){loc = i;};
		int getLoc(){return loc;};
		void setRSID(string r){rsid = r;};
		void setBPLOC(int l){bploc = l;};
		void setEnzyme(string e){enzyme = e;};
		void setProbeID(string p){probe_id = p;};
		void setChrom(string v){
			if(v == "X"){
				chrom = opts::_CHRX_;
			}
			else if(v == "Y"){
				chrom = opts::_CHRY_;
			}
			else if(v == "XY"){
				chrom = opts::_CHRXY_;
			}
			else{
				chrom = atoi(v.c_str());
			}
		};
		void setChrom(int v){chrom = v;};
		void setEnabled(bool v){enabled = v;};
		void setFlag(bool v){flag = v;};
		void setAllele1(string v){
       if (alleles.empty()) {
               alleles.resize(1, v);
            } else {
               alleles.at(0) = v;
            }
	  //  alleles.push_back(v);
		};//allele1 = v;};

		void setAllele2(string v){
       if (alleles.size() < 2) {
               alleles.resize(2, v);
            } else {
               alleles.at(1) = v;
            }
      //alleles.push_back(v);
		};//allele2 = v;};

    void setReferent(string ref){referent_allele = ref;}
    string getReferent(){return referent_allele;}
    int getReferentIndex(){return ref_all_index;}
    void setReferentIndex(int index){ref_all_index = index;}

		void addAllele(string v){
			alleles.push_back(v);
      // if(alleles.size() > 2){
      //  micro_sat.set(0);
      // }
		};

		bool isMicroSat(){
      return (alleles.size() > 2) ? true : false;
      // if(micro_sat[0]){
      //  return true;
      // }
      // return false;
		};

		int getAlleleLoc(string v){
			int i = -1;
			for(unsigned int j = 0; j < alleles.size(); j++){
				if(alleles[j] == v){
					i = j;
					return i;
				}
			}
			return i;
		};

		string getProbeID(){return probe_id;};
		string getEnzyme(){return enzyme;};
		string getRSID(){return rsid;};
		int getChrom(){return chrom;};
		int getBPLOC(){return bploc;};
		bool isEnabled(){return enabled;};
		bool isFlagged(){return flag;};
		string getAllele1(){
      // if(alleles.empty()) {
      //  this->setAllele1("");
      // }
      // return alleles.at(0);
			if(alleles.size() == 0) {
        return "";
			}
      return alleles[0];
		};//return allele1;};
		string getAllele2() {
      // if(alleles.size() < 2){
      //  this->setAllele2("");
      // }
      // return alleles.at(1);
			if(alleles.size() < 2){
        return "";
			}
      return alleles[1];
		};//return allele2;};

		int getNumAlleles(){
			return alleles.size();
		};
		vector<string> getAlleles(){
			return alleles;
		};
		string getAllele(int l){
			if(l >= (int)alleles.size()){
				cerr << l << " is > " << alleles.size() << endl;
			}
      // Possible solution for segfault.
      // return alleles.at(l);
      return alleles[l];
		};
		void resetAllele1(string a1){
			alleles[0] = a1;
		};
		void resetAllele2(string a2){
			alleles[1] = a2;
		};
		bool hasMAF(){
			return freqflag;
		};
		float getMAF(){
			return maf;
		};
		void setMAF(float v){
			maf = v;
		};
		void setFreqFlag(bool v){
			freqflag = v;
		};
		void assignDetail(string d, string v){
			details[d] = v;
		};
		string getDetailHeaders(){
			string value = "";
			map<string, string>::iterator iter;
			for(iter = details.begin(); iter != details.end(); iter++){
				value += iter->first + "\t";
			}
			if(value.length() > 0){
				value = value.erase(value.length() - 1, 1);
			}
			return value;
		};
		string getDetails(){
			string value = "";
			map<string, string>::iterator iter;
			for(iter = details.begin(); iter != details.end(); iter++){
				value += iter->second + "\t";
			}
			if(value.length() > 0){
				value = value.erase(value.length() - 1, 1);
				value = "\t" + value;
			}
			return value;
		};

};
};
#endif
