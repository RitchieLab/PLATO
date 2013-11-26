#ifndef CONCORDANCE_H
#define CONCORDANCE_H

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
#include "DataSet.h"
#include "StepOptions.h"

using namespace std;
namespace Methods{
class Concordance{
	static string stepname;
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
		string sample_error_file;
		string main_file;
		string error_file;
		string mismatch_file;

		vector<Sample*> check_samples;
		vector<Marker*> check_markers;
		vector<Family*> check_families;
		vector<int> check_marker_map;
		DataSet* check_data_set;

	public:
		Concordance(){
			families = NULL;
			markers = NULL;
			samples = NULL;
			check_data_set = new DataSet();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			_MARKERLIST_ = false;
			_STRATIFY_ = false;
			threshold = 0;
			order = 0;
			sample_error_file = "";
			main_file = "";
			error_file = "";
		};
		Concordance(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			samples = ds->get_samples();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			check_data_set = new DataSet();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			order = 0;
			sample_error_file = "";
			main_file = "";
			error_file = "";
		};
		Concordance(float thresh) : threshold(thresh){
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
			sample_error_file = "";
			main_file = "";
			error_file = "";
		};
		~Concordance(){
			check_samples.resize(0);
			check_families.resize(0);
			check_markers.resize(0);
			check_marker_map.clear();
			delete(check_data_set);
		};
		string get_sample_error_file(){return sample_error_file;}
		string get_main_file(){return main_file;}
		string get_error_file(){return error_file;}

		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
		};
		void FilterSummary();
		void setOptions(StepOptions o){
			options = o;
		};
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void calculate();
        void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
};

#endif
