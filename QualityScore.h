#ifndef QUALITYSCORE_H
#define QUALITYSCORE_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
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
#include "Marker.h"
#include "Sample.h"
#include "Family.h"
#include "Globals.h"
#include "Process.h"
using namespace std;


class QualityScore : public Process{
	private:
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;

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
		int order;

		float min_bin;
		float max_bin;
		float bin_space;
		map<float, int > qs_total_het;
		map<float, int > qs_total_maj;
		map<float, int > qs_total_min;
		vector<float> qs_ave;
		vector<int> qs_tot;

	public:
		QualityScore(){
			//families = NULL;
			//markers = NULL;
			markers = NULL;
			samples = NULL;
			families = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			threshold = 0;
			rank = 0;
			min_bin = 0.0f;
			max_bin = 0.5f;
			bin_space = 0.01f;
			order = 0;
		};
		QualityScore(float thresh) : threshold(thresh){
			markers = NULL;
			samples = NULL;
			families = NULL;
			marker_map = NULL;
			//families = NULL;
			//markers = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			min_bin = 0.0f;
			max_bin = 0.5f;
			bin_space = 0.01f;
			order = 0;
		};
		~QualityScore(){};
		//void process(Connection*, Families*, Markers*);
	//	void process(Families*, Markers*);
	//	void process(Connection*, vector<Sample*>*, vector<Family*>*, vector<Marker*>*){};
		void PrintSummary();
		void filter();
		void setThreshold(string s){threshold = std::atof(s.c_str());};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
		//void updateFamsMarks(Families* f, Markers* m){
		//	families = f;
	//		markers = m;
	//	};
		void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		int get_marker_loc(int);
				
};

#endif
