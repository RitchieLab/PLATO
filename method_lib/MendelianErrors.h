#ifndef MENDELIANERRORS_H
#define MENDELIANERRORS_H

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
#include "Family.h"
#include "Marker.h"
#include "Sample.h"
#include "DataSet.h"
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"

using namespace std;

namespace Methods{
class MendelianErrors{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		StepOptions options;

		int fam_thresh;
		int marker_thresh;
		void Tokenize(const string&, vector<string>&, const string&);
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
		string errors_file_name;
		string level2_file_name;

		vector<int> ferrors;
		vector< map<string,int> > fenzyme;
		vector<int> merrors;
		vector<int> serrors;
		vector< map<string,int> > senzyme;
		vector<vector<Marker*> > error_map;

		vector<Marker*> good_markers;


	public:
		MendelianErrors(){
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
		};
		MendelianErrors(DataSet* ds){
			data_set = ds;
			markers = ds->get_markers();
			families = ds->get_families();
			samples = ds->get_samples();
			marker_map = ds->get_marker_map();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			fam_thresh = -1;
			marker_thresh = -1;
			error_rate = 0.1;
			rank = 0;
			order = 0;
		};

		MendelianErrors(string thresh){
			vector<string> tokens;
			Tokenize(thresh, tokens, ":");
			if(tokens.size() == 2){
				fam_thresh = std::atoi(tokens.at(0).c_str());
				marker_thresh = std::atoi(tokens.at(1).c_str());
				cout << fam_thresh << ", " << marker_thresh << endl;
			}
			else{
				cerr << "Incorrect number of Mendelian Error thresholds.  Value should be #:#" << endl;
				exit(1);
			}
			fam_thresh = -1;
			marker_thresh = -1;
			markers = NULL;
			families = NULL;
			samples = NULL;
			marker_map = NULL;
			error_rate = 0.1;
			rank = 0;
			order = 0;
		};
		virtual ~MendelianErrors(){};

		string get_error_filename(){return errors_file_name;}
		string get_level2_filename(){return level2_file_name;}
		void PrintSummary();
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void setThreshold(string thresh);
		void setOptions(StepOptions o){
			options = o;
		};
		void filter_markers();
		void perform_evaluation(bool);
		void FilterSummary();
		void resetCounts();
		void calcThreshold();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
		void mid_process_output();
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void write_error(ofstream &out, Marker*, Sample*, Sample*, int);
		void zeroErrors();
		bool getPossibleGenos(Family*, Marker*);
		void determineZygotes(Family*, Marker*);
		int checkTypes(Family*, Marker*);
		void printError(ofstream &, Family*, Marker*);
		bool checkErrors(Family*, Marker*);
		bool find_match(vector<vector<int> >, vector<int>);
		void genotypeElimination(Family*, Marker*, int&);
		void genotypeElimination2(Family*, Marker*, int&);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void calculate(){perform_evaluation(false);};
		vector<int> getNumFamilyErrors(){return ferrors;};
		vector<int> getNumMarkerErrors(){return merrors;};
		vector<int> getNumSampleErrors(){return serrors;};
		vector<vector<Marker*> > getErrorMap(){return error_map;};
};
};
#endif
