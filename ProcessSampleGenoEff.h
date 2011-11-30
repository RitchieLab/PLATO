#ifndef PROCESSSAMPLEGENOEFF_H
#define PROCESSSAMPLEGENOEFF_H

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
#include <Marker.h>
#include <Family.h>
#include <Sample.h>
#include <Options.h>
#include <Globals.h>
#include "Process.h"
#include <StepOptions.h>
#include <SampleGenoEff.h>
#include <MethodException.h>
#include <DataSet.h>

using namespace std;

typedef vector<int> PERCENT;

class ProcessSampleGenoEff : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;
		//Markers* markers;
		//Families* families;
		float threshold;
		PERCENT per_cutoff;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		vector<int> zeros;
		vector<int> total;
		vector< map<string, int> > enzyme_zeros;
		vector< map<string, int> > enzyme_total;

	public:
		ProcessSampleGenoEff(){
			data_set = NULL;
			per_cutoff.push_back(100);
			per_cutoff.push_back(95);
			per_cutoff.push_back(90);
			per_cutoff.push_back(85);
			per_cutoff.push_back(80);
			per_cutoff.push_back(75);
			per_cutoff.push_back(70);
			per_cutoff.push_back(65);
			per_cutoff.push_back(60);
			per_cutoff.push_back(55);
			per_cutoff.push_back(50);
			families = NULL;
			markers = NULL;
			samples = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_samples = 0;
		};
		ProcessSampleGenoEff(float thresh) : threshold(thresh){
			data_set = NULL;
			per_cutoff.push_back(100);
			per_cutoff.push_back(95);
			per_cutoff.push_back(90);
			per_cutoff.push_back(85);
			per_cutoff.push_back(80);
			per_cutoff.push_back(75);
			per_cutoff.push_back(70);
			per_cutoff.push_back(65);
			per_cutoff.push_back(60);
			per_cutoff.push_back(55);
			per_cutoff.push_back(50);
			families = NULL;
			markers = NULL;
			samples = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_samples = 0;
		};
		~ProcessSampleGenoEff(){
		zeros.resize(0);
		total.resize(0);
		};
		//void process(Connection*, Families*, Markers*);
		//void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		PERCENT* getPerCutoff(){return &per_cutoff;};
		void process(DataSet*);
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumSamples(){return orig_num_samples;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
//        void updateFamsMarks(Families* f, Markers* m){
//			families = f;
//			markers = m;
//		};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};

#endif
