#ifndef SAMPLE_H
#define SAMPLE_H

#include <stdio.h>
#include <iostream>
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <list>
#include <bitset>
#include "Globals.h"
#include "General.h"
#include "Family.h"
#include "Trait.h"
#include "Covariate.h"

using namespace std;

namespace Methods{
class Family;
class Sample{
//case-control affected status: 0 = unaffected, > 0 = affected (inddx.status in the database)
	private:
		string famid;
		string id;
		int id_digit;
		string plate;
		string well;
		string dad;
		int dad_digit;
		string mom;
		int mom_digit;
		Sample* pdad;
		Sample* pmom;
		Sample* sibling;
		Sample* msibling;
		Sample* psibling;
		Sample* spouse;
		Family* pfam;
		bool sex; //true == male, false == female
		bool affected;
		bool enabled;
		bool excluded;
		bool flag;
		double phenotype;
		int pheno_loc;
		int loc;
		bool founder;
		int sysid;
		map<string,string> details;
		vector<Sample*> children;

        // 0/0 aone = true, atwo = true, amissing = true
		// A/A aone = false, atwo = false, amissing = false
		// T/T aone = true, atwo = true, amissing = false
		// A/T aone = false, atwo = true, amissing = false
		vector<bool> aone;
		vector<bool> atwo;
		vector<int> micro_sat_list;

		//vectors hold array location of microsat alleles in Marker.h
		//values of -1 -1 == zero call
		vector<int> abone;
		vector<int> abtwo;

		int tempabone;
		int tempabtwo;
		bool tempaone;
		bool tempatwo;

		//storage for deletion and homozygous steps
		vector<bool> roh;
		vector<char> deletion_cat;

		//possibly defunct; storage for level2 mendelian error checking
		map<int, vector<int> > micro_sat_alleles;
		vector<bool> aoneposs;
		vector<bool> atwoposs;
		vector<bool> amissingposs;
		vector<int> aboneposs;
		vector<int> abtwoposs;
		vector<bool> aonezyg;
		vector<bool> atwozyg;
		vector<int> abonezyg;
		vector<int> abtwozyg;
		vector<vector<int> > zygotes;


		bool isMEerror;
		vector<int> mesaved;

		vector<double> covs;
		vector<double> traits;
		vector<Trait*> trait_objs;
		vector<Covariate*> cov_objs;
		string grouping;

		///For QFAM model
		double T;
		double B;
		double W;

		int generation;

	public:
		Sample(){
			famid = id = plate = well = dad = mom = "";
			sibling = spouse = pdad = pmom = psibling = msibling = NULL;
			pfam = NULL;
			sex = affected = enabled = flag = false;
			phenotype = -1;
			aone.resize(0);
			atwo.resize(0);

			// add smd
			missing.resize(0);

			abone.resize(0);
			abtwo.resize(0);
			loc = 0;
			founder = false;
			isMEerror = false;
			id_digit = -1;
			dad_digit = -1;
			mom_digit = -1;
			excluded = false;
			T = B = W = 0;
			pheno_loc = -1;
			generation = 0;
		};
		~Sample(){};
        struct mysort{
		    bool operator()(const string s1, const string s2) const{
	            if(s1 == "X"){
		            return false;
	            }
	            else if(s2 == "X"){
	                return true;
	            }
	            else{
	                return atoi(s1.c_str()) < atoi(s2.c_str());
	            }
	        }
	    };

		int getGeneration(){return generation;}
		void setGeneration(int i){generation = i;}
		void incrGeneration(){generation++;}
		void decrGeneration(){generation--;}

		//For QFAM Model
		double getT(){return T;}
		void setT(double v){T = v;}
		double getW(){return W;}
		void setW(double v){W = v;}
		double getB(){return B;}
		void setB(double v){B = v;}

		void setSysid(int i ){sysid = i;}
		int getSysid(){return sysid;}

		string getFamIDOrig(){return famid;}
		string getIndOrig(){return id;}
		string getMomIDOrig(){return mom;}
		string getDadIDOrig(){return dad;}
		string toString();
		void setGrouping(string g){grouping = g;};
		string getGrouping(){return grouping;};
		string getFamID();
		void setFamID(string f){famid = f;};

		void assignDetail(string d, string v){
			details.at(d) = v;
		};
		void addCovariate(double v){covs.push_back(v);};
		void resizeCovariates(int c){covs.resize(c);};
		void setCovariate(double v, int l){covs.at(l) = v;};
		void set_covariate_vector(vector<double> v){covs = v;};
		double getCovariate(int i){return covs.at(i);};
		vector<double> getCovariateVector(){return covs;}
		void resizeTraits(int t){traits.resize(t);};
		void addTrait(double v){traits.push_back(v);};
		void setTrait(double v, int l){traits.at(l) = v;};
		void set_trait_vector(vector<double> v){traits = v;};
		double getTrait(int i){return traits.at(i);};
		vector<double> getTraitVector(){return traits;}

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

		void setIndDigit(int i){id_digit = i;};
		int getInd_digit(){
			if(id_digit == -1){
				id_digit = atoi(id.c_str());
			}
			return id_digit;
		};
		int getDadID_digit(){
			if(dad_digit == -1){
				dad_digit = atoi(dad.c_str());
			}
			return dad_digit;
		};
		int getMomID_digit(){
			if(mom_digit == -1){
				mom_digit = atoi(mom.c_str());
			}
			return mom_digit;
		};
		void setExcluded(bool b){excluded = b;};
		bool isExcluded(){return excluded;};
		void setParentDigits(){
			if(pdad != NULL){
				dad_digit = pdad->getInd_digit();
			}
			if(pmom != NULL){
				mom_digit = pmom->getInd_digit();
			}
		};
		void resizeDeletionCat(int s, char c){deletion_cat.resize(s, c);};
		char getDeletionCat(int s){return deletion_cat.at(s);};
		void setDeletionCat(int s, char c){deletion_cat.at(s) = c;};
		void resizeRoh(int s){roh.resize(s, false);};
		vector<bool> getRoh(){return roh;};
		bool getRohVal(int l){return roh.at(l);};
		void clearRoh(){roh.clear();};
		void setRoh(int l, bool v){roh.at(l) = v;};
		void setME2error(bool v){isMEerror = v;};
		bool IsMEerror(){return isMEerror;};
		void addMEsaved(int v){
			vector<int>::iterator found = find(mesaved.begin(), mesaved.end(), v);
			if(found == mesaved.end()){
				mesaved.push_back(v);
			}
		};
		vector<int> getMEsaved(){return mesaved;};
		int getMEsavedCount(){return mesaved.size();};
		void setMEsaved(vector<int> v){mesaved = v;};
		int getMEsaved(int l){return mesaved.at(l);};
		void clearMEsaved(){mesaved.clear();};
		string genoToString(bool ms, int loc, vector<string> alleles){
			string val = "";
			if(!ms){
				bool a1 = aone.at(loc);
				bool a2 = atwo.at(loc);
				bool a3 = missing.at(loc);
				if(a1 && a2 && a3){
					return "0/0";
				}
				if(!a1){
					val += alleles.at(0);
				}
				else{
					val += alleles.at(1);
				}
				if(!a2){
					val += "/" + alleles.at(0);
				}
				else{
					val += "/" + alleles.at(1);
				}
			}
			else{
				int a1 = micro_sat_alleles.at(loc).at(0);
				int a2 = micro_sat_alleles.at(loc).at(1);
				if(a1 == -1 && a2 == -1){
					return "0/0";
				}
				val += alleles.at(a1) + "/" + alleles.at(a2);
			}
			return val;

		};
		string possibleToString(vector<string> alleles){
			string val = "";
			for(unsigned int i = 0; i < aoneposs.size(); i++){
				if(aoneposs.at(i) && atwoposs.at(i) && amissingposs.at(i)){
					val += "0/0";
				}
				else{
					if(aoneposs.at(i) && atwoposs.at(i)){
						val += alleles.at(1) + "/" + alleles.at(1);
					}
					else if(!aoneposs.at(i) && !atwoposs.at(i)){
						val += alleles.at(0) + "/" + alleles.at(0);
					}
					else if(!aoneposs.at(i) && atwoposs.at(i)){
						val += alleles.at(0) + "/" + alleles.at(1);
					}
					else if(aoneposs.at(i) && !atwoposs.at(i)){
						val += alleles.at(1) + "/" + alleles.at(0);
					}
				}
				val += " ";
			}
			for(unsigned int i = 0; i < aboneposs.size(); i++){
				if(aboneposs.at(i) < (int)alleles.size() && abtwoposs.at(i) < (int)alleles.size()){
					val += alleles.at(aboneposs.at(i)) + "/" + alleles.at(abtwoposs.at(i)) + " ";
				}
				else{
					val += "(bad)" + getString<int>(aboneposs.at(i)) + "/" + getString<int>(abtwoposs.at(i)) + " ";
				}
			}
			return val;
		};

		string getDetails(){
			string value = "";
			map<string, string>::iterator iter;
			for(iter = details.begin(); iter != details.end(); iter++){
				value += iter->second + "\t";
			}
			if(value.length() > 0){
				value = value.erase(value.length() - 1, 1);
			}
			return value;
		};
		vector<vector<int> >* getZygotes(){return &zygotes;};
		void setZygotes(vector<vector<int> > v){zygotes = v;};
		void addZygotes(vector<vector<int> > v){
			for(unsigned int i = 0; i < v.size(); i++){
				vector<vector<int> >::iterator found = find(zygotes.begin(), zygotes.end(), v.at(i));
				if(found == zygotes.end()){
					zygotes.push_back(v.at(i));
				}
			}
		};
		void addMSPossible(int a1, int a2){
			aboneposs.push_back(a1);
			abtwoposs.push_back(a2);
		};
		void addPossible(bool a1, bool a2, bool am){
			aoneposs.push_back(a1);
			atwoposs.push_back(a2);
			amissingposs.push_back(am);
		};

		void remMSPossible(int loc){
			if(aboneposs.size() > 0){
				aboneposs.erase(aboneposs.begin() + loc);
				abtwoposs.erase(abtwoposs.begin() + loc);
			}
		};
		void remPossible(int loc){
			if(aoneposs.size() > 0){
				aoneposs.erase(aoneposs.begin() + loc);
				atwoposs.erase(atwoposs.begin() + loc);
				amissingposs.erase(amissingposs.begin() + loc);
			}
		};
		void clearPossible(){
			aoneposs.clear();
			atwoposs.clear();
			amissingposs.clear();
			aboneposs.clear();
			abtwoposs.clear();
		};
		vector<bool> getAonePossible(){
			return aoneposs;
		};
		vector<bool> getAtwoPossible(){
			return atwoposs;
		};
		vector<bool> getAmissingPossible(){
			return amissingposs;
		};
		vector<int> getAbonePossible(){
			return aboneposs;
		};
		vector<int> getAbtwoPossible(){
			return abtwoposs;
		};
		bool getAonePossibleVal(int l){return aoneposs.at(l);};
		bool getAtwoPossibleVal(int l){return atwoposs.at(l);};
		bool getAmissingPossibleVal(int l){return amissingposs.at(l);};
		int getAbonePossibleVal(int l){return aboneposs.at(l);};
		int getAbtwoPossibleVal(int l){return abtwoposs.at(l);};
		void setAonePossible(vector<bool> v){aoneposs = v;};
		void setAtwoPossible(vector<bool> v){atwoposs = v;};
		void setAmissingPossible(vector<bool> v){amissingposs = v;};
		void setAbonePossible(vector<int> v){aboneposs = v;};
		void setAbtwoPossible(vector<int> v){abtwoposs = v;};
		int validMSPossible(){
			int count = 0;
			for(unsigned int i = 0; i < aboneposs.size(); i++){
				if(aboneposs.at(i) != -1 && abtwoposs.at(i) != -1){
					count++;
				}
			}
			return count;
		};
		int validPossible(){
			int count = 0;
			for(unsigned int i = 0; i < aoneposs.size(); i++){
				if(!(aoneposs.at(i) && !atwoposs.at(i))){
					count++;
				}
			}
			return count;
		};

		void resizeAlleles(int c){
			aone.clear();
			aone.resize(c);
			atwo.clear();
			atwo.resize(c);

			// added smd
			missing.clear();
			missing.resize(c, false);
		};

		bool haveMicroSat(int b){
			for(unsigned int i = 0; i < micro_sat_list.size(); i++){
				if(micro_sat_list.at(i) == b){
					return true;
				}
			}
			return false;
		};
		int getMicroSat(int b){
			for(unsigned int i = 0; i < micro_sat_list.size(); i++){
				if(micro_sat_list.at(i) == b){
					return i;
				}
			}
			return -1;
		};
		void addMicroSat(int b){
			vector<int> a;
			a.resize(2);
			//micro_sat_alleles.at(b) = a;
			// changed to allow for setting of element
			// at throws exception if key doesn't exist in map	
			micro_sat_alleles[b]=a;
		};
		void addAbone(int m, int b){
			micro_sat_alleles.at(m).at(0) = b;
		};
		void addAbtwo(int m, int b){
			micro_sat_alleles.at(m).at(1) = b;
		};
		void addAbone(int b){
			abone.push_back(b);
		};
		void addAbtwo(int b){
			abtwo.push_back(b);
		};

		void addAone(bool b){
			aone.push_back(b);
		};
		void addAtwo(bool b){
			atwo.push_back(b);
		};
		void addAone(int i, bool b){
			aone.at(i) = b;
		};
		void addAtwo(int i, bool b){
			atwo.at(i) = b;
		};

		int getAoneSize(){
			return aone.size();
		};

		int getAtwoSize(){
			return atwo.size();
		};

		int getAboneSize(){
			return abone.size();
		};
		int getAbtwoSize(){
			return abtwo.size();
		};

		bool getAone(int i){
			return aone.at(i);
		};
		bool getAtwo(int i){
			return atwo.at(i);
		};

		int getAbone(int i){
			return micro_sat_alleles.at(i).at(0);
		};
		int getAbtwo(int i){
			return micro_sat_alleles.at(i).at(1);
		};

		void initializeTempAlleles(int loc){
			tempaone = aone.at(loc);
			tempatwo = atwo.at(loc);
		};

		void initializeTempBAlleles(int loc){
			tempabone = abone.at(loc);
			tempabtwo = abtwo.at(loc);
		};

		bool getTempAllele1(){
			return tempaone;
		};
		bool getTempAllele2(){
			return tempatwo;
		};

		int getTempBAllele1(){
			return tempabone;
		};
		int getTempBAllele2(){
			return tempabtwo;
		};

		void setTempAllele1(bool b){
			tempaone = b;
		};
		void setTempAllele2(bool b){
			tempatwo = b;
		};

		void setTempBAllele1(int i){
			tempabone = i;
		};
		void setTempBAllele2(int i){
			tempabtwo = i;
		};

		int getLoc(){return loc;};
		void setLoc(int l){loc = l;};

		bool isFounder(){return founder;};
		void setFounder(bool v){founder = v;};
		string getInd();
		void setInd(string v){id = v;};
		bool getSex(){return sex;};
		void setSex(bool v){sex = v;};
		string getSexAsString(){
			if(sex){
				return "M";
			}
			else{
				return "F";
			}
		};
		string getMomID();
		string getDadID();
		void setMomID(string v){mom = v;};
		void setDadID(string v){dad = v;};
		bool getAffected(){return affected;};
		void setAffected(bool v){affected = v;};
		Sample* getMom(){return pmom;};
		Sample* getDad(){return pdad;};
		void setMom(Sample* s){pmom = s;};
		void setDad(Sample* s){pdad = s;};
		void setSib(Sample* s){sibling = s;};
		void setMatSib(Sample* s){msibling = s;};
		void setPatSib(Sample* s){psibling = s;};
		Sample* getSib(){return sibling;};
		Sample* getMatSib(){return msibling;};
		Sample* getPatSib(){return psibling;};
		Family* getFamily(){return pfam;};
		void setFamily(Family* f){pfam = f;};
		vector<Sample*>* getChildren(){return &children;};
		void addChild(Sample* s){children.push_back(s);};
		string getPlate(){return plate;};
		void setPlate(string s){plate = s;};
		string getWell(){return well;};
		void setWell(string w){well = w;};
		bool isEnabled(){return enabled;};
		void setEnabled(bool val){enabled = val;};
		bool isFlagged(){return flag;};
		void setFlag(bool v){flag = v;};
		void setPheno(double v){phenotype = v;};
		void setUsePheno(int ploc){pheno_loc = ploc;}
		double getPheno(){
			if(pheno_loc >= 0 && (int)traits.size() > pheno_loc){
				return traits.at(pheno_loc);
			}
			return phenotype;
		};

		//returns trait.at(i) in trait list unless i >= traits size, then return phenotype
		double getPheno(int i){
			if(i >= 0 && (int)traits.size() > i){
				return traits.at(i);
			}
			return phenotype;
		}
		Sample* getLastChild();
template <class T>
 std::string getString(const T& t){
	        std::stringstream os;
			        os << t;
					        return os.str();
};

  // smd additions
  vector<bool> missing;

  unsigned int get_genotype(unsigned int marker_index){return aone.at(marker_index)+atwo.at(marker_index)+missing.at(marker_index);}
  void addAmissing(bool amiss){missing.push_back(amiss);}
  void addAmissing(int i, bool b){
    missing.at(i) = b;
  }
	bool getAmissing(int i){return missing.at(i);};
};
};

#endif
