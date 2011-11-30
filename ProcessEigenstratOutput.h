#ifndef PROCESSEIGENSTRATOUTPUT_H
#define PROCESSEIGENSTRATOUTPUT_H

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
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <EigenstratOutput.h>

using namespace std;
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessEigenstratOutput : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		//StepOptions options;

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

	public:
		ProcessEigenstratOutput(){
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
		ProcessEigenstratOutput(int thresh) : threshold(thresh){
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
		~ProcessEigenstratOutput(){};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = 0;
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
		//string map_allele(string);
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		int get_marker_loc(int);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void run(DataSetObject*){};
		void dump2db(){};
};
#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
#endif
