#ifndef GENDERCHECK_H
#define GENDERCHECK_H

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
#include "Marker.h"
#include "Family.h"
#include "Sample.h"
#include "Options.h"
#include "Globals.h"
#include "StepOptions.h"
#include "DataSet.h"
using namespace std;

namespace Methods{
class GenderCheck {
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
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		vector<Marker*> good_markers;

		vector<int> merrors;
		vector<int> shets;
		vector<int> mtotal;
		vector<int> stotal;
		vector< map<string, int> > senzyme_hets;
		vector< map<string, int> > senzyme_tot;

	public:
		GenderCheck(){
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
		};
		GenderCheck(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			orig_num_families = 0;
			orig_num_markers = 0;
			orig_num_samples = 0;
			ind_thresh = -1.0;
			marker_thresh = -1.0;
			error_rate = 0.1;
			rank = 0;
			order = 0;
		};

		GenderCheck(float thresh){
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
		};
		~GenderCheck(){};
		void setOptions(StepOptions o){
			options = o;
		};
		void PrintSummary();
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void calculate(){perform_evaluation(false);};
		void perform_evaluation(bool);
		void filter_markers();
		void setThreshold(string s);
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		void calcThreshold();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		void setOrder(int o){order = o;};

        void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		vector<int> getMarkerErrors(){return merrors;};
		vector<int> getSampleHets(){return shets;};
		vector<int> getMarkerTotalSamplesUsed(){return mtotal;};
		vector<int> getSampleTotalMarkersUsed(){return stotal;};
		vector< map<string, int> > getEnzymeHets(){return senzyme_hets;};
		vector< map<string, int> > getEnzymeTotals(){return senzyme_tot;};
		int getTotalSamplesUsedAtMarker(int i){
			if(i < (int)mtotal.size()){
				return mtotal[i];
			}
			else{
				return -1;
			}
		};
		int getNumErrorsAtMarker(int i){
			if(i < (int)merrors.size()){
				return merrors[i];
			}
			else{
				return -1;
			}
		};
		float getPercentMalesHetAtMarker(int i){
			if(i < (int)mtotal.size()){
				return (((float)merrors[i]/(float)mtotal[i]) * 100.0f);
			}
			else{
				return -1;
			}
		};
		int getNumHetsAtSample(int i){
			if(i < (int)shets.size()){
				return shets[i];
			}
			else{
				return -1;
			}
		};
		int getTotalMarkersUsedAtSample(int i){
			if(i < (int)stotal.size()){
				return stotal[i];
			}
			else{
				return -1;
			}
		};
		float getPercentHetAtSample(int i){
			if(i < (int)stotal.size()){
				return (((float)shets[i]/(float)stotal[i]) * 100.0f);
			}
			else{
				return -1;
			}
		};
};
};

#endif
