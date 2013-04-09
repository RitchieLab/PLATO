#ifndef PROCESSRUNTDT_H
#define PROCESSRUNTDT_H

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
#include <Sample.h>
//#include "TDTProcess.h"
#include <Options.h>
#include <StepOptions.h>
#include <RunTDT.h>
#include <DataSet.h>
#include <MethodException.h>

using namespace std;
using namespace Methods;

class ProcessRunTDT : public ProcessImpl<ProcessRunTDT>{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
	//	StepOptions options;
//		Markers* markers;
//		Families* families;
//		TDTProcess* TDT;
		int rank;
		float threshold;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_individuals;
		int order;

		vector<double> chi;
		vector<long double> pval;
		vector<int> fams_used;
		vector<float> maf;
		vector<double> trans;
		vector<double> untrans;

		//by group variables
		vector<vector<double> > gchi;
		vector<vector<long double> > gpval;
		vector<vector<int> > gfams_used;
		vector<vector<float> > gmaf;
		vector<vector<double> > gtrans;
		vector<vector<double> > guntrans;

		vector<Marker*> good_markers;

	public:
		ProcessRunTDT(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			//TDT = NULL;
			rank = 0;
			threshold = 0;
			orig_num_markers = 0;
			order = 0;
		};
		virtual ~ProcessRunTDT(){};
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s);
		void FilterSummary();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
		void setDBOUT(){_DBOUTPUT_ = true;};

		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
#endif
