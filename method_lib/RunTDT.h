#ifndef RUNTDT_H
#define RUNTDT_H

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
#include "Marker.h"
#include "Family.h"
#include "Globals.h"
#include "Sample.h"
#include "Options.h"
#include "DataSet.h"
#include "StepOptions.h"

using namespace std;

namespace Methods{
class RunTDT{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;
		int rank;
		float threshold;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int orig_num_markers;
		int order;

		vector<double> chi;
		vector<long double> pval;
		vector<int> fams_used;
		vector<float> maf;
		vector<double> trans;
		vector<double> untrans;

		//method
		double chi_one;
		long double pval_one;
		int fams_used_one;
		float maf_one;
		double trans_one;
		double untrans_one;
		double odds_ratio_one;
		double ci_l_one;
		double ci_u_one;

	public:
		RunTDT(){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			rank = 0;
			threshold = 0;
			orig_num_markers = 0;
			order = 0;
		};

		RunTDT(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			markers = ds->get_markers();
			samples = ds->get_samples();
			marker_map = ds->get_marker_map();
			rank = 0;
			threshold = 0;
			orig_num_markers = 0;
			order = 0;
		};

		~RunTDT(){};
		void PrintSummary();
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void setThreshold(string s);
		void setOptions(StepOptions o){
			options = o;
		};
		void FilterSummary();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};

		void calculate(Marker*);
		void calculate(int m){calculate((*markers).at(m));};
		double getChi(){return chi_one;};
		long double getPval(){return pval_one;};
		int getFamsUsed(){return fams_used_one;};
		double getTransmitted(){return trans_one;};
		double getUntransmitted(){return untrans_one;};
		double getOddsRatio(){return odds_ratio_one;};
		double getConfIntervalLower(){return ci_l_one;};
		double getConfIntervalUpper(){return ci_u_one;};
		void docalc(Sample* samp, Sample* dad, Sample* mom, Marker* mark, double& trans1, double &trans2, int &fams);
		void calculate(Marker* mark, vector<Sample*> gsamps);



};
};
#endif
