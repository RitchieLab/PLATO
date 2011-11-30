#ifndef LD_H
#define LD_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
//#include "Markers.h"
//#include "Families.h"
#include "Globals.h"
#include "Options.h"
#include "General.h"
#include "StepOptions.h"
#include "DataSet.h"

using namespace std;


namespace Methods{
class LD{// : public Process{
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;
		//Markers* markers;
		//Families* families;
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
		int run_chr;
		int run_start;
		int run_end;


	public:
		LD(){
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
			run_chr = 0;
			run_start = -1;
			run_end = -1;
		};
		LD(float thresh) : threshold(thresh){
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
			run_chr = 0;
			run_start = -1;
			run_end = -1;
		};
		~LD(){
		};
		//void process(Connection*, Families*, Markers*);
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		void setOptions(StepOptions o){
			options = o;
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		//void updateFamsMarks(Families* f, Markers* m){
		//	families = f;
		//	markers = m;
		//};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		//void process(Families*, Markers*);
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
        void calculate(DataSet* ds){
			data_set = ds;
			process(ds->get_samples(), ds->get_families(), ds->get_markers(), ds->get_marker_map());
		};
		void setOrder(int o){order = o;};
		void calcSetMeanVariance(vector<Marker*>, vector<double>&, vector<vector<double> > &);
		//vector< vector<double> > svd_inverse(vector<vector<double> > &);
		//void svdcmp(vector<vector<double> > &, vector<double> &, vector<vector<double> > &);
		//double pythag(const double, const double);
		//double SQR(double);
		vector<bool> vif_prune(vector<vector<double> >, double);
		vector<int> getChromosomeMarkerRange(vector<Marker*>*, int);
		void setMarkerRange();
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void calcLDStatistics();
		vector<double> correlation2SNP(Marker*, Marker*);
};
};
#endif

