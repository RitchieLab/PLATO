#ifndef HOMOZYGOUS_H
#define HOMOZYGOUS_H

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
//#include "Process.h"
#include "Options.h"
#include "General.h"
#include "StepOptions.h"
#include "DataSet.h"
#include "MethodException.h"

using namespace std;

namespace Methods{
class Homozygous {//: public Process{
	static string stepname;
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
		vector<int> homoaffcount;
		vector<int> homounaffcount;
		vector<int> homoallcount;
		vector<int> homominallcount;
		vector<int> homomajallcount;

		vector<Marker*> good_markers;

		vector<string> filenames;


	public:
		Homozygous(){
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
		};

		Homozygous(DataSet* ds){
			data_set = ds;
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
		};

		Homozygous(float thresh) : threshold(thresh){
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
		};
		~Homozygous(){
		};
		//void process(Connection*, Families*, Markers*);
		//
		vector<string> get_filenames(){return filenames;}
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		void setOptions(StepOptions o){
			options = o;
		};
		void set_parameters(StepOptions* o){
			options = *o;
		};
		void resetDataSet(DataSet* ds){
			data_set = ds;
		};
		void calculate(){process(data_set->get_samples(), data_set->get_families(), data_set->get_markers(), data_set->get_marker_map());};
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
        void setOrder(int o){order = o;};
		double checkSequence(int, int, int&, int&, double&);
//		vector<int> findCommonRoh(Sample*, vector<Marker*>*, vector<int>*);
		void process_wgha();
		vector<Marker*> enabledMarkers();
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void perform_homozyg_permutations();
		Sample* findRandomSample(map<Sample*, bool>);
		vector<int> getHomoUnaffCount(){return homounaffcount;};
		vector<int> getHomoAffCount(){return homoaffcount;};
		vector<int> getHomoMajAllCount(){return homomajallcount;};
		vector<int> getHomoMinAllCount(){return homominallcount;};
		vector<int> getHomoAllCount(){return homoallcount;};
};
};
#endif
