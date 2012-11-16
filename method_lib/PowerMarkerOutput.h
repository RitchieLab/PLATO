#ifndef POWERMARKEROUTPUT_H
#define POWERMARKEROUTPUT_H

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
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"
#include "DataSet.h"

using namespace std;


namespace Methods{
class PowerMarkerOutput{// : public Process{
	private:
		DataSet* data_set;
		static string stepname;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;

//		Markers* markers;
//		Families* families;
		int threshold;	
//		PERCENT per_cutoff;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_individuals;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;
		vector<string> filenames;

	public:
		PowerMarkerOutput(){
//			per_cutoff.push_back(100);
//			per_cutoff.push_back(95);
//			per_cutoff.push_back(90);
//			per_cutoff.push_back(85);
//			per_cutoff.push_back(80);
//			per_cutoff.push_back(75);
//			per_cutoff.push_back(70);
//			per_cutoff.push_back(65);
//			per_cutoff.push_back(60);
//			per_cutoff.push_back(55);
//			per_cutoff.push_back(50);
			data_set = NULL;
			families = NULL;
			markers = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		PowerMarkerOutput(int thresh) : threshold(thresh){
//			per_cutoff.push_back(100);
//			per_cutoff.push_back(95);
//			per_cutoff.push_back(90);
//			per_cutoff.push_back(85);
//			per_cutoff.push_back(80);
//			per_cutoff.push_back(75);
//			per_cutoff.push_back(70);
//			per_cutoff.push_back(65);
//			per_cutoff.push_back(60);
//			per_cutoff.push_back(55);
//			per_cutoff.push_back(50);
			data_set = NULL;
			families = NULL;
			markers = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		~PowerMarkerOutput(){};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		vector<string> get_filenames(){return filenames;}
		void filter();
		void process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm);
		void setThreshold(string s){
			options.setUp(s);
			if(options.haveCenterCodes()){
				options.readCenterFile(options.getCenterFile());
			}
			//	threshold = 0;
		};
		void setOptions(StepOptions o){
			options = o;
		};
		void calculate(DataSet* ds){
			data_set = ds;
			process(ds->get_samples(), ds->get_families(), ds->get_markers(), ds->get_marker_map());
		};
//		PERCENT* getPerCutoff(){return &per_cutoff;};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		Sample* find_sample(string, string);
		bool find_marker(string);
        void setOrder(int o){order = o;};
//		void updateFamsMarks(Families* f, Markers* m){
//			families = f;
//			markers = m;
//		};
		void setDBOUT(){_DBOUTPUT_ = true;};
		int map_sex(char);
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		int get_marker_loc(int);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
};
#endif
