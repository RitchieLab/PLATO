#ifndef PROCESSCACONCHISQ_H
#define PROCESSCACONCHISQ_H

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
//#include <ChiSquare.h>
#include <Globals.h>
#include <General.h>
#include "Process.h"
#include <CaConChisq.h>
#include <MultComparison.h>
#include <Options.h>
#include <StepOptions.h>
#include <MethodException.h>
#include <DataSet.h>
using namespace std;
using namespace Methods;

class ProcessCaConChisq : public ProcessImpl<ProcessCaConChisq>{
	static string stepname;
	private:
		DataSet* data_set;
		//StepOptions options;

//		Markers* markers;
//		Families* families;
		float threshold;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		vector<double> chi_geno;
		vector<double> chi_allele;
		vector<double> chi_arm;
		vector<double> pval_arm;
		vector<double> pval_allele;
		vector<double> pval_geno;
		vector<double> pval_geno_exact;
		vector<double> pval_allele_exact;
		vector<double> odds_ratio;
		vector<double> ci_l;
		vector<double> ci_u;
		vector<int> geno_df;
		vector<int> allele_df;
		vector<int> arm_df;

		//groups
		double gchi_geno;
		double gchi_allele;
		double gchi_arm;
		double gpval_arm;
		double gpval_allele;
		double gpval_geno;
		double gpval_geno_exact;
		double gpval_allele_exact;
		double godds_ratio;

		vector<Marker*> good_markers;

	public:
		ProcessCaConChisq(){
			data_set = NULL;
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
		};
		ProcessCaConChisq(float thresh) : threshold(thresh){
			data_set = NULL;
			rank = 0;
			order =0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
		};
#ifdef PLATOLIB
		ProcessCaConChisq(string, int, Database*);
#endif

		virtual ~ProcessCaConChisq();

		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getRank(){return rank;};
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void FilterSummary();
		void resize(int);
		void setRank(int r){rank = r;};
		void setOrder(int o){order = o;};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void setThreshold(string s){options.setUp(s);};
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
