#ifndef PROCESSMITOCHECK_H
#define PROCESSMITOCHECK_H

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
#include <Options.h>
#include <Globals.h>
#include "Process.h"
#include <StepOptions.h>
#include <MitoCheck.h>
#include <MethodException.h>
#include <DataSet.h>

using namespace std;

class ProcessMitoCheck : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;
		float ind_thresh;
		float marker_thresh;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;
		float error_rate;
		void Tokenize(const string&, vector<string>&, const string&);
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		vector<int> merrors;
		vector<int> serrors;
		vector<vector<int> > error_map;

	public:
		ProcessMitoCheck(){
			data_set = NULL;
			orig_num_families = 0;
			orig_num_markers = 0;
			orig_num_samples = 0;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			ind_thresh = -1.0;
			marker_thresh = -1.0;
			error_rate = 0.1;
			rank = 0;
			order = 0;
			overwrite = true;
		};
		ProcessMitoCheck(float thresh){
			data_set = NULL;
			orig_num_families = 0;
			orig_num_markers = 0;
			orig_num_samples = 0;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			ind_thresh = -1.0;
			marker_thresh = -1.0;
			error_rate = 0.1;
			rank = 0;
			order =0;
			overwrite = true;
		};
		~ProcessMitoCheck(){};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
//		void perform_evaluation(Connection*, bool);
		void filter_markers();
		void setThreshold(string s);
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		void setOrder(int o){order = o;};
//        void updateFamsMarks(Families* f, Markers* m){
//		    families = f;
//		    markers = m;
//		};

        void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};

#endif
