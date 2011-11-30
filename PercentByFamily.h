#ifndef PERCENTBYFAMILY_H
#define PERCENTBYFAMILY_H

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
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <Sample.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <MethodException.h>

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
using namespace PlatoLib;
#endif

//typedef vector<int> PERCENT;

class PercentByFamily : public Process{
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
		//PERCENT per_cutoff;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		vector<int> fzeros;
		vector<int> ftotal;
		vector<int> mzeros;
		vector<int> mtotal;
		vector< map<string, int> > enzyme_zeros;
		vector< map<string, int> > enzyme_total;

	public:
		PercentByFamily(){
			/*per_cutoff.push_back(100);
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
			*/
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			rank = 0;
			order = 0;
		};
		PercentByFamily(float thresh) : threshold(thresh){
			/*per_cutoff.push_back(100);
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
			*/
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			rank = 0;
			order =0;
		};
		~PercentByFamily(){};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s){threshold = std::atof(s.c_str());};
		//PERCENT* getPerCutoff(){return &per_cutoff;};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_samples;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
//		void updateFamsMarks(Families* f, Markers* m){
//			families = f;
//			markers = m;
//		};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
        void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void run(DataSetObject*){};
		void dump2db(){};
};

#endif
