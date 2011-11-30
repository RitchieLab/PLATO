#ifndef HWEQUILIBRIUM_H
#define HWEQUILIBRIUM_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "AlleleFrequency.h"
#include "ChiSquare.h"
#include "Marker.h"
#include "Family.h"
#include "Sample.h"
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"

using namespace std;
namespace Methods{
class HWEquilibrium{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;
		float threshold;	
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;	
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool useoverall;
		bool overwrite;
		int order;

		AlleleFrequency* af;
		float hw_O;
		float hw_OM;
		float hw_OF;
		float hw_P;
		float hw_PM;
		float hw_PF;
		float hw_C;
		float hw_CM;
		float hw_CF;
		float hw_Ca;
		float hw_CaF;
		float hw_CaM;
		float hw_Con;
		float hw_ConM;
		float hw_ConF;

		//HWEPT
		double geno_chi;
		double geno_pval;
		double allele_chi;
		double allele_p;
		double ca_pval;
		double con_pval;


		vector<float> hwO;
		vector<float> hwP;
		vector<float> hwPM;
		vector<float> hwPD;
		vector<float> hwC;
		vector<float> hwCM;
		vector<float> hwCF;
		vector<float> hwCa;
		vector<float> hwCaF;
		vector<float> hwCaM;
		vector<float> hwCon;
		vector<float> hwConF;
		vector<float> hwConM;
			
		
	public:
		HWEquilibrium(){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			rank = 0;
			_MARKERLIST_ = false;
			_DBOUTPUT_ = false;
			_STRATIFY_ = false;
			order = 0;
			orig_num_markers = 0;
			useoverall = false;
		};
		HWEquilibrium(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			samples = ds->get_samples();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			af = new AlleleFrequency(ds);
			order = 0;
			rank = 0;
			orig_num_markers = 0;
			useoverall = false;
		};
		HWEquilibrium(float thresh) : threshold(thresh){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			rank = 0;
			_MARKERLIST_ = false;
			_DBOUTPUT_ = false;
			_STRATIFY_ = false;
			order = 0;
			orig_num_markers = 0;
			useoverall = false;
		};
		~HWEquilibrium(){
		};
		void resetDataSet(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			samples = ds->get_samples();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			af = new AlleleFrequency(ds);
			af->setOptions(options);
			orig_num_markers = 0;
		};
		void PrintSummary();
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void setThreshold(string s){
			options.setUp(s);
		};
		void setOptions(StepOptions o){
			options = o;
			if(options.doRandomChild() || options.doAll() || options.doAllChildren()){
				useoverall = true;
			}
			af->setOptions(options);
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};

		float calcHW(float, float, int, int, int, int);
		float calcHW_exact(int, int, int);
		void resize(int i);
		void doFilter(Marker*);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void processHWEPT();
		void calculateHWEPT(Marker*);
		vector<double> hwePT(vector<bool>, Marker*);
		void doFilterHWEPT(Marker*);

		void calculate(int m){calculate((*markers).at(m));};
		void calculate(Marker* m);
		
		float getOverall(){return hw_O;};
		float getOverallMale(){return hw_OM;};
		float getOverallFemale(){return hw_OF;};
		float getParental(){return hw_P;};
		float getParentalMale(){return hw_PM;};
		float getParentalFemale(){return hw_PF;};
		float getChildren(){return hw_C;};
		float getChildrenMale(){return hw_CM;};
		float getChildrenFemale(){return hw_CF;};
		float getCase(){return hw_Ca;};
		float getCaseFemale(){return hw_CaF;};
		float getCaseMale(){return hw_CaM;};
		float getControl(){return hw_Con;};
		float getControlMale(){return hw_ConM;};
		float getControlFemale(){return hw_ConF;};
		//HWEPT
		double getGenotypicChi(){return geno_chi;};
		double getGenotypicPval(){return geno_pval;};
		double getAllelicChi(){return allele_chi;};
		double getAllelicPval(){return allele_p;};
		double getCasePvalHWEPT(){return ca_pval;};
		double getControlPvalHWEPT(){return con_pval;};
		float getCaseGeneralHWEPT(){return hw_Ca;};
		float getControlGeneralHWEPT(){return hw_Con;};

};
};

#endif
