#ifndef PROCESSBEAGLEOUTPUT_H
#define PROCESSBEAGLEOUTPUT_H

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
#include <Marker.h>
#include <Sample.h>
#include <Family.h>
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <MethodException.h>
#include <BEAGLEOutput.h>

using namespace std;
using namespace Methods;
//define the PlatoLib namespace for use with Plato as a library
#ifdef PLATOLIB
namespace PlatoLib
{
#endif
class ProcessBEAGLEOutput : public ProcessImpl<ProcessBEAGLEOutput>{
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
		int order;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		string projectPath;

	public:
		ProcessBEAGLEOutput(){
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
		ProcessBEAGLEOutput(int thresh) : threshold(thresh){
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
#ifdef PLATOLIB
		ProcessBEAGLEOutput(string, int, Database*, string);
#endif
		virtual ~ProcessBEAGLEOutput(){};
//		void process(Families*, Markers*);
//		void process(Connection*, Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s){
			options.setUp(s);
			//threshold = 0;
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
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
        bool hasIncExc(){return options.doIncExcludedSamples();};
		#ifdef PLATOLIB
			void run(DataSetObject*);
			void dump2db();
			void create_tables();
		#endif
};
#ifdef PLATOLIB
};
#endif
#endif
