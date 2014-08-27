#ifndef CMH_H
#define CMH_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "Globals.h"
#include "Options.h"
#include "General.h"
#include "StepOptions.h"
#include "DataSet.h"

using namespace std;

namespace Methods{
class CMH{
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;
		float threshold;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_individuals;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;
		void calcCMH(Marker*);
		void calcMantelHaenszel_2x2xK(Marker*);
		vector<double> calcMantelHaenszel_IxJxK(vector<int> & X, vector<int> & Y,
				vector<int> & Z);
		vector<double> calcMantelHaenszel_ORD(vector<int> & X, vector<int> & Y,
				vector<int> & Z);

		double chisq;
		double chisq_bd;
		double pval;
		double pval_bd;
		double OR_upper;
		double OR_lower;
		double OR;
		vector<bool> samp_flags;

	public:
		CMH(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			_MARKERLIST_ = false;
			_STRATIFY_ = false;
			threshold = 0;
			order = 0;
			chisq = 0;
			pval = 0;
			chisq_bd = 0;
			pval_bd = 0;
			OR_upper = 0;
			OR_lower = 0;
			OR = 0;
		};
		CMH(float thresh) : threshold(thresh){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			_MARKERLIST_ = false;
			_STRATIFY_ = false;
			order = 0;
			chisq = 0;
			pval = 0;
			chisq_bd = 0;
			pval_bd = 0;
			OR_upper = 0;
			OR_lower = 0;
			OR = 0;
		};
		~CMH(){
		};
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
		};
		void setOptions(StepOptions o){
			options = o;
			if(!options.doClusterFile()){
				throw MethodException("-cluster-file is a required argument for CMH.");
			}
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
        void calculate(Marker*);
        void calculate(int i){calculate(data_set->get_locus(i));};
        void resetDataSet(DataSet* ds);
		void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void flagSamples();
		double get_chisq(){return chisq;}
		double get_pval(){return pval;}
		double getOR(){return OR;}
		double getOR_upper(){return OR_upper;}
		double getOR_lower(){return OR_lower;}
		double get_chisqbd(){return chisq_bd;}
		double get_pvalbd(){return pval_bd;}
};
};

#endif
