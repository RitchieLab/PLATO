#ifndef CACONCHISQ_H
#define CACONCHISQ_H

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
#include "AlleleFrequency.h"
#include "ChiSquare.h"
#include "Marker.h"
#include "Family.h"
#include "Sample.h"
#include "Globals.h"
//#include "Process.h"
#include "Options.h"
#include "StepOptions.h"
#include "DataSet.h"
using namespace std;

namespace Methods{
class CaConChisq{// : public Process{
	static string stepname;
	private:
		DataSet* data_set;

		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;

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

		AlleleFrequency* af;
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

		//method
		double chi_geno_one;
		double chi_allele_one;
		double chi_arm_one;
		double pval_arm_one;
		double pval_allele_one;
		double pval_geno_one;
		double pval_geno_exact_one;
		double pval_allele_exact_one;
		double odds_ratio_one;
		double ci_l_one;
		double ci_u_one;
		int allele_exact_df;
		int geno_df;
		int allele_df;
		int arm_df;
		map<string, double> gchi_geno_one;
		map<string, double> gchi_allele_one;
		map<string, double> gchi_arm_one;
		map<string, double> gpval_arm_one;
		map<string, double> gpval_allele_one;
		map<string, double> gpval_geno_one;
		map<string, double> gpval_geno_exact_one;
		map<string, double> gpval_allele_exact_one;
		map<string, double> godds_ratio_one;
		map<string, double> gci_l_one;
		map<string, double> gci_u_one;

	public:
		CaConChisq(){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			rank = 0;
			order = 0;
		orig_num_markers = 0;
		orig_num_families = 0;
		orig_num_samples = 0;
		};

		CaConChisq(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			markers = ds->get_markers();
			samples = ds->get_samples();
			marker_map = ds->get_marker_map();
			rank = 0;
			af = new AlleleFrequency(ds);
			af->setRank(rank);
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
		};

		CaConChisq(float thresh) : threshold(thresh){
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			rank = 0;
			order =0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_samples = 0;
		};
		~CaConChisq(){
			//if(af != NULL){
			//delete(af);
			//}
		};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
			//threshold = std::atof(s.c_str());
		};
		void setOptions(StepOptions o){
			options = o;
			if(options.doGroupFile()){
				options.readGroups(samples);
			}
			af->setOptions(options);
		};

		void calculate(Marker* mark);
		void calculate(int m){calculate((*markers)[m]);};
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
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

		double calcChiGeno(int, int, int, int, int, int);
	//	double calcChiGeno_exact(int, int, int);
		//double calcChiAllele(vector<vector<int> >);
		//double calcChiAlleleFisher(int, int, int, int);

		void resize(int i);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};

		//method
		double getGenotypicChi(){return chi_geno_one;};
		int getGenotypicDF(){return geno_df;}
		double getAllelicChi(){return chi_allele_one;};
		int getAllelicDF(){return allele_df;}
		double getArmitageChi(){return chi_arm_one;};
		int getArmitageDF(){return arm_df;}
		double getArmitagePval(){return pval_arm_one;};
		double getAllelicPval(){return pval_allele_one;};
		double getGenotypicPval(){return pval_geno_one;};
		double getGenotypicExactPval(){return pval_geno_exact_one;};
		double getAllelicExactPval(){return pval_allele_exact_one;};
		double getOddsRatio(){return odds_ratio_one;};
		double getConfIntervalLower(){return ci_l_one;};
		double getConfIntervalUpper(){return ci_u_one;};
		map<string, double> getGroupGenotypicChi(){return gchi_geno_one;};
		double getGroupGenotypicChi(string v){return gchi_geno_one[v];};
		map<string, double> getGroupAllelicChi(){return gchi_allele_one;};
		double getGroupAllelicChi(string v){return gchi_allele_one[v];};
		map<string, double> getGroupArmitageChi(){return gchi_arm_one;};
		double getGroupArmitageChi(string v){return gchi_arm_one[v];};
		map<string, double> getGroupArmitagePval(){return gpval_arm_one;};
		double getGroupArmitagePval(string v){return gpval_arm_one[v];};
		map<string, double> getGroupAllelicPval(){return gpval_allele_one;};
		double getGroupAllelicPval(string v){return gpval_allele_one[v];};
		map<string, double> getGroupGenotypicPval(){return gpval_geno_one;};
		double getGroupGenotypicPval(string v){return gpval_geno_one[v];};
		map<string, double> getGroupGenotypicExactPval(){return gpval_geno_exact_one;};
		double getGroupGenotypicExactPval(string v){return gpval_geno_exact_one[v];};
		map<string, double> getGroupAllelicExactPval(){return gpval_allele_exact_one;};
		double getGroupAllelicExactPval(string v){return gpval_allele_exact_one[v];};
		map<string, double> getGroupOddsRatio(){return godds_ratio_one;};
		double getGroupOddsRatio(string v){return godds_ratio_one[v];};
		map<string, double> getGroupConfIntervalLower(){return gci_l_one;};
		double getGroupConfIntervalLower(string v){return gci_l_one[v];};
		map<string, double> getGroupConfIntervalUpper(){return gci_u_one;};
		double getGroupConfIntervalUpper(string v){return gci_u_one[v];};


    // reset data
    void resetData(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			markers = ds->get_markers();
			samples = ds->get_samples();
			marker_map = ds->get_marker_map();
			rank = 0;
			af = new AlleleFrequency(ds);
			af->setRank(rank);
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
	  }


};
};

#endif
