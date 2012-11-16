#ifndef ALLELEINFO_H
#define ALLELEINFO_H

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <occi.h>
#include <string>
#include <list>
#include <map>
#include "Globals.h"
#include "ChiSquare.h"

using namespace oracle::occi;
using namespace std;

class AlleleInfo{
	private:
		string allele1;
		string allele2;
		float case_major_freq_f;
		float case_minor_freq_f;
		float case_major_freq_m;
		float case_minor_freq_m;
		float control_major_freq_f;
		float control_minor_freq_f;
		float control_major_freq_m;
		float control_minor_freq_m;
		
        float parent_major_freq;
        float parent_minor_freq;
        float father_major_freq;
        float mother_major_freq;
        float father_minor_freq;
        float mother_minor_freq;
        float child_major_freq;
        float child_minor_freq;
		float child_major_freq_f;
		float child_minor_freq_f;
		float child_major_freq_m;
		float child_minor_freq_m;
        int parent_major_count;
        int parent_minor_count;
        int father_major_count;
        int mother_major_count;
        int father_minor_count;
        int mother_minor_count;
        int child_major_count;
        int child_minor_count;
        int child_major_count_f;
        int child_minor_count_f;
        int child_major_count_m;
        int child_minor_count_m;

		int case_major_count_f;
		int case_minor_count_f;
		int case_major_count_m;
		int case_minor_count_m;
		int control_major_count_f;
		int control_minor_count_f;
		int control_major_count_m;
		int control_minor_count_m;
		
        int parent_het_count;
        int father_het_count;
        int mother_het_count;
        int child_het_count;
        int child_het_count_f;
        int child_het_count_m;

		int case_het_count_f;
		int case_het_count_m;
		int control_het_count_f;
		int control_het_count_m;
		
        int parent_homo1_count;
        int parent_homo2_count;
        int father_homo1_count;
        int mother_homo1_count;
        int father_homo2_count;
        int mother_homo2_count;
        int child_homo1_count;
        int child_homo2_count;
        int child_homo1_count_f;
        int child_homo2_count_f;
        int child_homo1_count_m;
        int child_homo2_count_m;
		
		int case_homo1_count_f;
		int case_homo2_count_f;
		int case_homo1_count_m;
		int case_homo2_count_m;
		int control_homo1_count_f;
		int control_homo2_count_f;
		int control_homo1_count_m;
		int control_homo2_count_m;
		
        float parent_hw;
        float father_hw;
        float mother_hw;
        float child_hw;
        float child_hw_f;
        float child_hw_m;

		float case_hw_f;
		float case_hw_m;
		float control_hw_m;
		float control_hw_f;
		
        int parents_used;
        int fathers_used;
        int mothers_used;
        int children_used;
        int children_used_f;
        int children_used_m;

		int cases_used_f;
		int cases_used_m;
		int controls_used_f;
		int controls_used_m;
		
	string major;
        void initialize(){
            parent_major_freq = 0.0;
            parent_minor_freq = 0.0;
            father_major_freq = 0.0;
            mother_major_freq = 0.0;
            father_minor_freq = 0.0;
            mother_minor_freq = 0.0;
            child_major_freq = 0.0;
            child_minor_freq = 0.0;
            child_major_freq_f = 0.0;
            child_minor_freq_f = 0.0;
            child_major_freq_m = 0.0;
            child_minor_freq_m = 0.0;

			case_major_freq_f = 0.0;
			case_major_freq_m = 0.0;
			case_minor_freq_f = 0.0;
			case_minor_freq_m = 0.0;
			control_major_freq_f = 0.0;
			control_major_freq_m = 0.0;
			control_minor_freq_f = 0.0;
			control_minor_freq_m = 0.0;

            parent_major_count = 0;
            parent_minor_count = 0;
            father_major_count = 0;
            mother_major_count = 0;
            father_minor_count = 0;
            mother_minor_count = 0;
            child_major_count = 0;
            child_minor_count = 0;
            child_major_count_f = 0;
            child_minor_count_f = 0;
            child_major_count_m = 0;
            child_minor_count_m = 0;

			case_major_count_f = 0;
			case_major_count_m = 0;
			case_minor_count_f = 0;
			case_minor_count_m = 0;
			control_major_count_f = 0;
			control_major_count_m = 0;
			control_minor_count_f = 0;
			control_minor_count_m = 0;
			
            parent_het_count = 0;
            father_het_count = 0;
            mother_het_count = 0;
            child_het_count = 0;
            child_het_count_f = 0;
            child_het_count_m = 0;
			
			case_het_count_f = 0;
			case_het_count_m = 0;
			control_het_count_f = 0;
			control_het_count_m = 0;
			
            parent_homo1_count = 0;
            father_homo1_count = 0;
            mother_homo1_count = 0;
            father_homo2_count = 0;
		    mother_homo2_count = 0;
            parent_homo2_count = 0;
            child_homo1_count = 0;
            child_homo2_count = 0;
            child_homo1_count_f = 0;
            child_homo2_count_f = 0;
            child_homo1_count_m = 0;
            child_homo2_count_m = 0;

			case_homo1_count_f = 0;
			case_homo1_count_m = 0;
			case_homo2_count_f = 0;
			case_homo2_count_m = 0;
			control_homo1_count_f = 0;
			control_homo1_count_m = 0;
			control_homo2_count_f = 0;
			control_homo2_count_m = 0;
			
            parent_hw = 0.0;
            father_hw = 0.0;
            mother_hw = 0.0;
            child_hw = 0.0;
            child_hw_f = 0.0;
            child_hw_m = 0.0;

			case_hw_f = 0.0;
			case_hw_m = 0.0;
			control_hw_f = 0.0;
			control_hw_m = 0.0;
			
            parents_used = 0;
            fathers_used = 0;
            mothers_used = 0;
            children_used = 0;
            children_used_f = 0;
            children_used_m = 0;

			cases_used_f = 0;
			cases_used_m = 0;
			controls_used_m = 0;
			controls_used_f = 0;
        };
		
			
	public:
		AlleleInfo(){
			initialize();
		};
		AlleleInfo(string a1, string a2){
			allele1 = a1;
			allele2 = a2;
			major = a1;
			initialize();
		};
		
		~AlleleInfo(){};
		void merge(allele_struct_data);
		void process(string, string, string);
		void process(bool, bool, string);
		void setMajorAllele(){
			if(parent_minor_count > parent_major_count){
				major = allele2;
			}
			else{
				major = allele1;
			}
		}
		string getMajorAllele(){return major;};
		float getParentMajorFreq(){
			if(parent_minor_freq > parent_major_freq){
				return parent_minor_freq;
			}
			return parent_major_freq;
		};
		float getParentMajorFreqAct(){
			return parent_major_freq;};
        float getFatherMajorFreq(){
			if(father_minor_freq > father_major_freq){
				return father_minor_freq;
			}
			return father_major_freq;
		};
        float getFatherMajorFreqAct(){
			return father_major_freq;};
        float getMotherMajorFreq(){
			if(mother_minor_freq > mother_major_freq){
				return mother_minor_freq;
			}
			return mother_major_freq;};
        float getMotherMajorFreqAct(){
			return mother_major_freq;};
		float getParentMinorFreq(){
			if(parent_minor_freq > parent_major_freq){
				return parent_major_freq;
			}
			return parent_minor_freq;};
		float getParentMinorFreqAct(){
			return parent_minor_freq;};
        float getFatherMinorFreq(){
			if(father_minor_freq > father_major_freq){
				return father_major_freq;
			}
			return father_minor_freq;};
        float getFatherMinorFreqAct(){
			return father_minor_freq;};
        float getMotherMinorFreq(){
			if(mother_minor_freq > mother_major_freq){
				return mother_major_freq;
			}
			return mother_minor_freq;};
        float getMotherMinorFreqAct(){
			return mother_minor_freq;};
		float getChildMajorFreq(){
			if(child_minor_freq > child_major_freq){
				return child_minor_freq;
			}
			return child_major_freq;};
		float getChildMajorFreqAct(){
			return child_major_freq;};
		float getChildMinorFreq(){
			if(child_minor_freq > child_major_freq){
				return child_major_freq;
			}
			return child_minor_freq;};
		float getChildMinorFreqAct(){
			return child_minor_freq;};
		float getChildMajorFreq_F(){
			if(child_minor_freq_f > child_major_freq_f){
				return child_minor_freq_f;
			}
			return child_major_freq_f;};
		float getChildMajorFreqAct_F(){
			return child_major_freq_f;};
		float getChildMinorFreq_F(){
			if(child_minor_freq_f > child_major_freq_f){
				return child_major_freq_f;
			}
			return child_minor_freq_f;};
		float getChildMinorFreqAct_F(){
			return child_minor_freq_f;};
		float getChildMajorFreq_M(){
			if(child_minor_freq_m > child_major_freq_m){
				return child_minor_freq_m;
			}
			return child_major_freq_m;};
		float getChildMajorFreqAct_M(){
			return child_major_freq_m;};
		float getChildMinorFreq_M(){
			if(child_minor_freq_m > child_major_freq_m){
				return child_major_freq_m;
			}
			return child_minor_freq_m;};
		float getChildMinorFreqAct_M(){
			return child_minor_freq_m;};

		float getCaseMajorFreq_F(){
			if(case_minor_freq_f > case_major_freq_f){
				return case_minor_freq_f;
			}
			return case_major_freq_f;};
		float getCaseMajorFreqAct_F(){
			return case_major_freq_f;};
		float getCaseMinorFreq_F(){
			if(case_minor_freq_f > case_major_freq_f){
				return case_major_freq_f;
			}
			return case_minor_freq_f;};
		float getCaseMinorFreqAct_F(){
			return case_minor_freq_f;};
		
		float getCaseMajorFreq_M(){
			if(case_minor_freq_m > case_major_freq_m){
				return case_minor_freq_m;
			}
			return case_major_freq_m;};
		float getCaseMajorFreqAct_M(){
			return case_major_freq_m;};
		float getCaseMinorFreq_M(){
			if(case_minor_freq_m > case_major_freq_m){
				return case_major_freq_m;
			}
			return case_minor_freq_m;};
		float getCaseMinorFreqAct_M(){
			return case_minor_freq_m;};

		float getControlMajorFreq_M(){
			if(control_minor_freq_m > control_major_freq_m){
				return control_minor_freq_m;
			}
			return control_major_freq_m;};
		float getControlMajorFreqAct_M(){
			return control_major_freq_m;};
		float getControlMinorFreq_M(){
			if(control_minor_freq_m > control_major_freq_m){
				return control_major_freq_m;
			}
			return control_minor_freq_m;};
		float getControlMinorFreqAct_M(){
			return control_minor_freq_m;};
		
		float getControlMajorFreq_F(){
			if(control_minor_freq_f > control_major_freq_f){
				return control_minor_freq_f;
			}
			return control_major_freq_f;};
		float getControlMajorFreqAct_F(){
			return control_major_freq_f;};
		float getControlMinorFreq_F(){
			if(control_minor_freq_f > control_major_freq_f){
				return control_major_freq_f;
			}
			return control_minor_freq_f;};
		float getControlMinorFreqAct_F(){
			return control_minor_freq_f;};


		
		int getParentMajorCount(){
			if(parent_minor_count > parent_major_count){
				return parent_minor_count;
			}
			return parent_major_count;};
		int getParentMajorCountAct(){
			return parent_major_count;};
        int getFatherMajorCount(){
			if(father_minor_count > father_major_count){
				return father_minor_count;
			}
			return father_major_count;};
        int getFatherMajorCountAct(){
			return father_major_count;};
        int getMotherMajorCount(){
			if(mother_minor_count > mother_major_count){
				return mother_minor_count;
			}
			return mother_major_count;};
        int getMotherMajorCountAct(){
			return mother_major_count;};
		int getParentMinorCount(){
			if(parent_minor_count > parent_major_count){
				return parent_major_count;
			}
			return parent_minor_count;};
		int getParentMinorCountAct(){
			return parent_minor_count;};
        int getFatherMinorCount(){
			if(father_minor_count > father_major_count){
				return father_major_count;
			}
			return father_minor_count;};
        int getFatherMinorCountAct(){
			return father_minor_count;};
        int getMotherMinorCount(){
			if(mother_minor_count > mother_major_count){
				return mother_major_count;
			}
			return mother_minor_count;};
        int getMotherMinorCountAct(){
			return mother_minor_count;};
		int getChildMajorCount(){
			if(child_minor_count > child_major_count){
				return child_minor_count;
			}
			return child_major_count;};
		int getChildMajorCountAct(){
			return child_major_count;};
		int getChildMinorCount(){
			if(child_minor_count > child_major_count){
				return child_major_count;
			}
			return child_minor_count;};
		int getChildMinorCountAct(){
			return child_minor_count;};
		int getChildMajorCount_F(){
			if(child_minor_count_f > child_major_count_f){
				return child_minor_count_f;
			}
			return child_major_count_f;};
		int getChildMajorCountAct_F(){
			return child_major_count_f;};
		int getChildMinorCount_F(){
			if(child_minor_count_f > child_major_count_f){
				return child_major_count_f;
			}
			return child_minor_count_f;};
		int getChildMinorCountAct_F(){
			return child_minor_count_f;};
		int getChildMajorCount_M(){
			if(child_minor_count_m > child_major_count_m){
				return child_minor_count_m;
			}
			return child_major_count_m;};
		int getChildMajorCountAct_M(){
			return child_major_count_m;};
		int getChildMinorCount_M(){
			if(child_minor_count_m > child_major_count_m){
				return child_major_count_m;
			}
			return child_minor_count_m;};
		int getChildMinorCountAct_M(){
			return child_minor_count_m;};

		//case/control
		int getCaseMajorCount_F(){
			if(case_minor_count_f > case_major_count_f){
				return case_minor_count_f;
			}
			return case_major_count_f;};
		int getCaseMajorCountAct_F(){
			return case_major_count_f;};
		int getCaseMinorCount_F(){
			if(case_minor_count_f > case_major_count_f){
				return case_major_count_f;
			}
			return case_minor_count_f;};
		int getCaseMinorCountAct_F(){
			return case_minor_count_f;};
		int getCaseMajorCount_M(){
			if(case_minor_count_m > case_major_count_m){
				return case_minor_count_m;
			}
			return case_major_count_m;};
		int getCaseMajorCountAct_M(){
			return case_major_count_m;};
		int getCaseMinorCount_M(){
			if(case_minor_count_m > case_major_count_m){
				return case_major_count_m;
			}
			return case_minor_count_m;};
		int getCaseMinorCountAct_M(){
			return case_minor_count_m;};
		int getControlMajorCount_F(){
			if(control_minor_count_f > control_major_count_f){
				return control_minor_count_f;
			}
			return control_major_count_f;};
		int getControlMajorCountAct_F(){
			return control_major_count_f;};
		int getControlMinorCount_F(){
			if(control_minor_count_f > control_major_count_f){
				return control_major_count_f;
			}
			return control_minor_count_f;};
		int getControlMinorCountAct_F(){
			return control_minor_count_f;};
		int getControlMajorCount_M(){
			if(control_minor_count_m > control_major_count_m){
				return control_minor_count_m;
			}
			return control_major_count_m;};
		int getControlMajorCountAct_M(){
			return control_major_count_m;};
		int getControlMinorCount_M(){
			if(control_minor_count_m > control_major_count_m){
				return control_major_count_m;
			}
			return control_minor_count_m;};
		int getControlMinorCountAct_M(){
			return control_minor_count_m;};
		
		
		int getParentHetCount(){return parent_het_count;};
        int getFatherHetCount(){return father_het_count;};
        int getMotherHetCount(){return mother_het_count;};
		int getChildHetCount(){return child_het_count;};
		int getChildHetCount_F(){return child_het_count_f;};
		int getChildHetCount_M(){return child_het_count_m;};
		//casecontrol
		int getCaseHetCount_F(){return case_het_count_f;};
		int getCaseHetCount_M(){return case_het_count_m;};
		int getControlHetCount_F(){return control_het_count_f;};
		int getControlHetCount_M(){return control_het_count_m;};
		
		int getParentHomo1Count(){return parent_homo1_count;};
        int getFatherHomo1Count(){return father_homo1_count;};
        int getMotherHomo1Count(){return mother_homo1_count;};
		int getParentHomo2Count(){return parent_homo2_count;};
        int getFatherHomo2Count(){return father_homo2_count;};
        int getMotherHomo2Count(){return mother_homo2_count;};
		int getChildHomo1Count(){return child_homo1_count;};
		int getChildHomo2Count(){return child_homo2_count;};
		int getChildHomo1Count_F(){return child_homo1_count_f;};
		int getChildHomo2Count_F(){return child_homo2_count_f;};
		int getChildHomo1Count_M(){return child_homo1_count_m;};
		int getChildHomo2Count_M(){return child_homo2_count_m;};
		//casecontrol
		int getCaseHomo1Count_F(){return case_homo1_count_f;};
		int getCaseHomo2Count_F(){return case_homo2_count_f;};
		int getCaseHomo1Count_M(){return case_homo1_count_m;};
		int getCaseHomo2Count_M(){return case_homo2_count_m;};
		int getControlHomo1Count_F(){return control_homo1_count_f;};
		int getControlHomo2Count_F(){return control_homo2_count_f;};
		int getControlHomo1Count_M(){return control_homo1_count_m;};
		int getControlHomo2Count_M(){return control_homo2_count_m;};
		
		void calcFreqs();
		float getABSParentMinorFreq();
        float getABSFatherMinorFreq();
        float getABSMotherMinorFreq();
		void reset(){
            parent_major_freq = parent_minor_freq =
	            case_major_freq_f = case_major_freq_m = case_minor_freq_f = case_minor_freq_m = 
	            control_major_freq_f = control_major_freq_m = control_minor_freq_f = control_minor_freq_m = 
	            child_major_freq = child_minor_freq = child_major_freq_f = child_major_freq_m = child_minor_freq_f = child_minor_freq_m = 
	            father_major_freq = father_minor_freq =
	            mother_major_freq = mother_minor_freq = 0.0;
            parent_major_count = parent_minor_count =
	            father_major_count = father_minor_count =
		        mother_major_count = mother_minor_count =
		        case_major_count_f = case_major_count_m = case_minor_count_f = case_minor_count_m = 
		        control_major_count_f = control_major_count_m = control_minor_count_f = control_minor_count_m = 
		        child_major_count = child_minor_count = child_major_count_f = child_major_count_m = child_minor_count_f = child_minor_count_m = 
		        case_het_count_f = case_het_count_m = 
		        control_het_count_f = control_het_count_m = 
		        parent_het_count = child_het_count = child_het_count_f = child_het_count_m = 
		        parent_homo1_count = parent_homo2_count =
		        case_homo1_count_f = case_homo1_count_m = case_homo2_count_f = case_homo2_count_m = 
		        control_homo1_count_f = control_homo1_count_m = control_homo2_count_f = control_homo2_count_m = 
		        child_homo1_count = child_homo2_count = child_homo1_count_f = child_homo1_count_m = child_homo2_count_f = child_homo2_count_m = 0;
		};
		float calcHW(float,float,int,int,int,int);
		void setParentHW(float val){parent_hw = val;};
		void setChildHW(float val){child_hw = val;};
		void setChildHW_F(float val){child_hw_f = val;};
		void setChildHW_M(float val){child_hw_m = val;};
        //casecontrol
		void setCaseHW_F(float val){case_hw_f = val;};
		void setCaseHW_M(float val){case_hw_m = val;};
		void setControlHW_F(float val){control_hw_f = val;};
		void setControlHW_M(float val){control_hw_m = val;};
		
		void setFatherHW(float val){father_hw = val;};
        void setMotherHW(float val){mother_hw = val;};
		float getParentHW(){return parent_hw;};
        float getFatherHW(){return father_hw;};
        float getMotherHW(){return mother_hw;};
		float getChildHW(){return child_hw;};
		float getChildHW_F(){return child_hw_f;};
		float getChildHW_M(){return child_hw_m;};
		//casecontrol
		float getCaseHW_F(){return case_hw_f;};
		float getCaseHW_M(){return case_hw_m;};
		float getControlHW_F(){return control_hw_f;};
		float getControlHW_M(){return control_hw_m;};
		
		float _subchisqrprob(int,float);
		float _subuprob(float);
        void setAlleles(string a1, string a2){allele1 = a1; allele2 = a2;};
		float getParentHetExp();
        float getFatherHetExp();
        float getMotherHetExp();
		float getChildHetExp();
		float getChildHetExp_F();
		float getChildHetExp_M();
        //casecontrol
		float getCaseHetExp_F();
		float getControlHetExp_F();
		float getCaseHetExp_M();
		float getControlHetExp_M();
		
		float getParentHomo1Exp();
        float getParentHomo2Exp();
        float getFatherHomo1Exp();
        float getFatherHomo2Exp();
        float getMotherHomo1Exp();
        float getMotherHomo2Exp();
        float getChildHomo1Exp();
        float getChildHomo2Exp();
        float getChildHomo1Exp_F();
        float getChildHomo2Exp_F();
        float getChildHomo1Exp_M();
        float getChildHomo2Exp_M();
        //casecontrol
        float getCaseHomo1Exp_F();
        float getCaseHomo2Exp_F();
        float getCaseHomo1Exp_M();
        float getCaseHomo2Exp_M();
        float getControlHomo1Exp_F();
        float getControlHomo2Exp_F();
        float getControlHomo1Exp_M();
        float getControlHomo2Exp_M();
		
		int getParentsUsed(){return parents_used;};
        int getFathersUsed(){return fathers_used;};
        int getMothersUsed(){return mothers_used;};
        int getChildrenUsed(){return children_used;};
        int getChildrenUsed_F(){return children_used_f;};
        int getChildrenUsed_M(){return children_used_m;};
		//casecontrol
        int getCasesUsed_F(){return cases_used_f;};
        int getCasesUsed_M(){return cases_used_m;};
        int getControlsUsed_F(){return controls_used_f;};
        int getControlsUsed_M(){return controls_used_m;};

		string getAllele1(){return allele1;};
		string getAllele2(){return allele2;};
		float calcHW_exact(int, int, int);
};

#endif
