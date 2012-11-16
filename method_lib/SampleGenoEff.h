#ifndef SAMPLEGENOEFF_H
#define SAMPLEGENOEFF_H

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
#include "Sample.h"
#include "Options.h"
#include "Globals.h"
//#include "Process.h"
#include "DataSet.h"
#include "StepOptions.h"
using namespace std;

namespace Methods{
	typedef vector<int> PERCENT;

class SampleGenoEff{
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

		int zeros_one;
		int total_one;

	
	public:
		SampleGenoEff(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
		};
		SampleGenoEff(){
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
		SampleGenoEff(float thresh) : threshold(thresh){
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
		~SampleGenoEff(){
		zeros.resize(0);
		total.resize(0);
		};
		//void process(Connection*, Families*, Markers*);
		//void process(Families*, Markers*);
		void resetDataSet(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
		};

		void PrintSummary();
		void filter();
		void filterOne(int s);
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		void setOptions(StepOptions o){
			options = o;
		};
		PERCENT* getPerCutoff(){return &per_cutoff;};
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void calculate(int s){calculate((*samples)[s]);};
		void calculate(Sample*);
		int getZeros(){return zeros_one;};
		int getTotal(){return total_one;};
		double getPercent(){return (1 - ((double)zeros_one/(double)total_one)) * 100.0f;};
			
//			(1 - (double)((double)zeros_one/(double)total_one)*100.0f);};
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
};
#endif
