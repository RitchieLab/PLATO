#ifndef PROCESSMENDELIANERRORS_H
#define PROCESSMENDELIANERRORS_H

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
#include "Process.h"
#include <Family.h>
#include <Marker.h>
#include <Sample.h>
#include "Chrom.h"
#include <Globals.h>
#include <Options.h>
#include <StepOptions.h>
#include <MendelianErrors.h>
#include <MethodException.h>
#include <DataSet.h>

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessMendelianErrors : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		//StepOptions options;

//		Markers* markers;
//		Families* families;
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
		string defaultinsert;
		string sampleinsert;
		string pedigreeinsert;
		string error_filename;
		string level2_filename;
		string projectPath;

		vector<int> ferrors;
		vector< map<string,int> > fenzyme;
		vector<int> merrors;
		vector<int> serrors;
		vector< map<string,int> > senzyme;
		vector<vector<Marker*> > error_map;
		vector<Marker*> good_markers;

	public:
		ProcessMendelianErrors() : Process(){
			data_set = NULL;
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
		ProcessMendelianErrors(string thresh){
			data_set = NULL;
			vector<string> tokens;
			Tokenize(thresh, tokens, ":");
			if(tokens.size() == 2){
				fam_thresh = std::atoi(tokens[0].c_str());
				marker_thresh = std::atoi(tokens[1].c_str());
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
		ProcessMendelianErrors(string, int, Database*, string);
		virtual ~ProcessMendelianErrors(){};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string thresh);
		void filter_markers();
		void FilterSummary();
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void zeroErrors();
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
