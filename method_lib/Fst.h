#ifndef FST_H
#define FST_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "AlleleFrequency.h"
#include "HWEquilibrium.h"
#include "Globals.h"
#include "Options.h"
#include "General.h"
#include "DataSet.h"
#include "StepOptions.h"

using namespace std;

namespace Methods{
class Fst{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;
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

		AlleleFrequency af;
		HWEquilibrium hwe;

		double fst;
		double fstrh;
		double fsthm;

		vector<int> zeros;
		vector<int> total;
		vector<int> casezeros;
		vector<int> casetotal;
		vector<int> controlzeros;
		vector<int> controltotal;
		map<string, vector<int> > groupzeros;
		map<string, vector<int> > grouptotal;

		int zeros_one;
		int total_one;
		int casezeros_one;
		int casetotal_one;
		int controlzeros_one;
		int controltotal_one;
		map<string, int> groupzeros_one;
		map<string, int> grouptotal_one;

	public:
		Fst(){
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
			fst = 0;
			fstrh = 0;
			fsthm = 0;
		};
		Fst(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			_MARKERLIST_ = false;
			_STRATIFY_ = false;
			threshold = 0;
			order = 0;
			if(options.doGroupFile() && samples != NULL){
				options.readGroups(samples);
			}
			af.resetDataSet(data_set);
			af.setOptions(options);


			fst = 0;
			fstrh = 0;
			fsthm = 0;
		};
		Fst(float thresh) : threshold(thresh){
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
			fst = 0;
			fstrh = 0;
			fsthm = 0;
		};
		~Fst(){
		zeros.resize(0);
		total.resize(0);
		casezeros.resize(0);
		casetotal.resize(0);
		controlzeros.resize(0);
		controltotal.resize(0);
		};
		void resetDataSet(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			if(options.doGroupFile() && samples != NULL){
				options.readGroups(samples);
			}
			af.resetDataSet(data_set);
			af.setOptions(options);

		};

		double getFst(){return fst;}
		double getFstRH(){return fstrh;}
		double getFstHM(){return fsthm;}
		void PrintSummary();
		void filter();
		void filterOne(int);
		void set_parameters(StepOptions* opts){
			options = *opts;
			if(options.doGroupFile() && samples != NULL){
				options.readGroups(samples);
			}
			if(options.getDoGroupFreq()){
				options.readGroupFrequencies();
			}
		};
		void setThreshold(string s){
			options.setUp(s);
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
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
        void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void calcOne(int);
		void calculate(int m){calcOne(m);};
		int getZeros(){return zeros_one;};
		int getTotal(){return total_one;};
		double getPercent(){ return (double)((1 - ((double)zeros_one/(double)total_one))*100.0f);};
		int getCaseZeros(){return casezeros_one;};
		int getCaseTotal(){return casetotal_one;};
		double getCasePercent(){ return (double)((1 - ((double)casezeros_one/(double)casetotal_one))*100.0f);};
		int getControlZeros(){return controlzeros_one;};
		int getControlTotal(){return controltotal_one;};
		double getControlPercent(){ return (double)((1 - ((double)controlzeros_one/(double)controltotal_one))*100.0f);};
		map<string, int> getGroupZeros(){return groupzeros_one;};
		map<string, int> getGroupTotal(){return grouptotal_one;};
};
};
#endif
