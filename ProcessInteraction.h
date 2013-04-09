#ifndef PROCESSINTERACTION_H
#define PROCESSINTERACTION_H

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
#include <Interactions.h>


using namespace std;
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessInteraction : public ProcessImpl<ProcessInteraction>{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		
		int threshold;
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
		ProcessInteraction(){
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

		ProcessInteraction(int thresh) : threshold(thresh){
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
		virtual ~ProcessInteraction(){};
		
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = 0;
		};
		
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		Sample* find_sample(string, string);
		bool find_marker(string);
    void setOrder(int o){order = o;};
    void setDBOUT(){_DBOUTPUT_ = true;};
    
    void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		int get_marker_loc(int);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		#ifdef PLATOLIB
			void run(DataSetObject*);
			void dump2db();
			void create_tables();
		#endif    
};
#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
#endif
