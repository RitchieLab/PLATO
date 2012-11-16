#ifndef PROCESSLOGREG_H
#define PROCESSLOGREG_H

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
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <General.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <LogisticRegression.h>
#include <cdflib.h>

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessLogReg : public Process{
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		//StepOptions options;
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
		int run_chr;
		int run_start;
		int run_end;
		string defaultinsert;

	public:
		ProcessLogReg(){
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
			run_chr = 0;
			run_start = -1;
			run_end = -1;
		};
		ProcessLogReg(float thresh) : threshold(thresh){
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
			run_chr = 0;
			run_start = -1;
			run_end = -1;
		};
#ifdef PLATOLIB
		ProcessLogReg(string, int, Database*);
#endif
		~ProcessLogReg(){};
		//void process(Connection*, Families*, Markers*);
		void PrintSummary();
		void filter();
		void doFilter(Methods::Marker*, double);
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		//void process(Families*, Markers*);
		void process(DataSet*);
        void setOrder(int o){order = o;};
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
