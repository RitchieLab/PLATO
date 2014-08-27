#ifndef KINSHIP_H
#define KINSHIP_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
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

typedef struct data{
	int father;
	int mother;
	int generation;
}parents;

struct info{
	map<int, parents> ped;
	int Nsize;
};

class Kinship{
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

		map<string, double> coefficients;


		vector<int> keep;
		vector< vector<int> > list;
		int i,j,F,fam;
		map<int, parents>::iterator iter, iterend;
	public:
		Kinship(){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			//TDT = NULL;
			rank = 0;
			threshold = 0;
			orig_num_markers = 0;
			order = 0;
		};

		Kinship(DataSet* ds){
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

		~Kinship(){};
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
		void resetDataSet(DataSet* ds){data_set = ds;}
		void calculate(Family*);
		void calculate(int m){calculate(data_set->get_pedigree(m));};
		double getChi(){return chi_one;};
		long double getPval(){return pval_one;};
		int getFamsUsed(){return fams_used_one;};
		double getTransmitted(){return trans_one;};
		double getUntransmitted(){return untrans_one;};
		double getOddsRatio(){return odds_ratio_one;};
		double getConfIntervalLower(){return ci_l_one;};
		double getConfIntervalUpper(){return ci_u_one;};


		void check(Sample*, Sample*, Family*);
		double phi2(Sample*, Sample*, Family*);
		void create_generation(Family*);
		int make_generation(Sample*, Sample*, Sample*, Family*);
		map<string, double> getCoefficients(){return coefficients;}

};
};
#endif
