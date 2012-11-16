#ifndef SUPERLINKOUTPUT_H
#define SUPERLINKOUTPUT_H

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
//#include "Markers.h"
//#include "Families.h"
#include "Marker.h"
#include "Sample.h"
#include "Family.h"
#include "Globals.h"
#include "Options.h"
#include "AlleleFrequency.h"
#include "StepOptions.h"
#include "DataSet.h"

using namespace std;


namespace Methods{
class SuperlinkOutput{// : public Process{
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
		SuperlinkOutput(){
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
			samples = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order =0;
		};
		SuperlinkOutput(int thresh) : threshold(thresh){
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
			samples = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		~SuperlinkOutput(){};
//		void process(Families*, Markers*);
//		void process(Connection*, Families*, Markers*);
		void PrintSummary();
		vector<string> get_filenames(){return filenames;}
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void setThreshold(string s){
			options.setUp(s);
			if(options.havePenetrance()){
				options.readPenetranceFile(options.getPenetranceFile());
			}
			else{
				opts::printLog(stepname + " requires a penetrance file to be included using the -penetrance-file option!\n");
				throw MethodException(stepname + " requires a penetrance file to be included using the -penetrance-file option!\n");
			}
			//threshold = 0;
		};
		void setOptions(StepOptions o){
			options = o;
			if(options.havePenetrance()){
				options.readPenetranceFile(options.getPenetranceFile());
			}
			else{
				opts::printLog(stepname + " requires a penetrance file to be included using the -penetrance-file option!\n");
				throw MethodException(stepname + " requires a penetrance file to be included using the -penetrance-file option!\n");
			}
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
        void setOrder(int o){order = o;};
//		void updateFamsMarks(Families* f, Markers* m){
//			families = f;
//			markers = m;
//		};
		void setDBOUT(){_DBOUTPUT_ = true;};
		int map_sex(char);
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
};
#endif
