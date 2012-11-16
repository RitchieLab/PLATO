#ifndef DELETIONS_H
#define DELETIONS_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <string>
#include <list>
//#include "Process.h"
#include "Family.h"
#include "Marker.h"
#include "Sample.h"
//#include "Chrom.h"
#include "Globals.h"
#include "Options.h"
#include "DataSet.h"
#include "StepOptions.h"

using namespace std;
namespace Methods{

class Deletions{// : public Process{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;

//		Markers* markers;
//		Families* families;
		int fam_thresh;
		int marker_thresh;
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

		string default_filename;
		string density_filename;

		vector<int> ferrors;
		vector< map<string,int> > fenzyme;
		vector<int> merrors;
		vector<int> serrors;
		vector< map<string,int> > senzyme;
		vector<vector<int> > error_map;
		
	public:
		Deletions(){// : Process(){
			markers = NULL; 
			samples = NULL;
			families = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			fam_thresh = -1;
			marker_thresh = -1;
			error_rate = 0.1;	
			rank = 0;
			order = 0;
			default_filename = "";
			density_filename = "";
		};
		Deletions(DataSet* ds){
			data_set = ds;
			markers = ds->get_markers();
			families = ds->get_families();
			marker_map = ds->get_marker_map();
			samples = ds->get_samples();
			rank = 0;
			order =0;
			default_filename = "";
			density_filename = "";
		};
		Deletions(string thresh){
			vector<string> tokens;
			fam_thresh = -1;
			marker_thresh = -1;
			markers = NULL;
			families = NULL;
			samples = NULL;
			marker_map = NULL;
			error_rate = 0.1;	
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			default_filename = "";
			density_filename = "";
		};
		virtual ~Deletions(){};
		string get_default_filename(){return default_filename;}
		string get_density_filename(){return density_filename;}
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void setOptions(StepOptions o){
			options = o;
		};
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void calculate(){perform_evaluation(false);};
		void setThreshold(string thresh);
		void filter_markers();
		void perform_evaluation(bool);
//		void perform_evaluation(Connection*, bool);
		void FilterSummary();
		void calcThreshold();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
//		void updateFamsMarks(Families* f, Markers* m){
//		    families = f;
//		    markers = m;
//		};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void zeroErrors(); 
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		};
		};
#endif
