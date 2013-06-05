#ifndef ALLELEFREQUENCY_H
#define ALLELEFREQUENCY_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <map>
#include "Marker.h"
#include "Family.h"
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"
#include "DataSet.h"
using namespace std;

namespace Methods{
class AlleleFrequency{
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
		bool useoverall;
		bool overwrite;
		int order;

		vector<bool> sample_flags;

		//overall groups
		map<string, int> ga1_count, ga2_count, ga1_homo_count, ga2_homo_count, ga12_count;

		//overall (no criteria...raw counts)
		int a1_count, a2_count, a1_homo_count, a2_homo_count, a12_count;
		//overall male
		int a1_countM, a2_countM, a1_homo_countM, a2_homo_countM, a12_countM;
		//overall female
		int a1_countF, a2_countF, a1_homo_countF, a2_homo_countF, a12_countF;
		//founder
		int a1_countP, a2_countP, a1_homo_countP, a2_homo_countP, a12_countP;
		//founder male
		int a1_countPM, a2_countPM, a1_homo_countPM, a2_homo_countPM, a12_countPM;
		//founder female
		int a1_countPF, a2_countPF, a1_homo_countPF, a2_homo_countPF, a12_countPF;
		//child
		int a1_countC, a2_countC, a1_homo_countC, a2_homo_countC, a12_countC;
		//child male
		int a1_countCM, a2_countCM, a1_homo_countCM, a2_homo_countCM, a12_countCM;
		//child female
		int a1_countCF, a2_countCF, a1_homo_countCF, a2_homo_countCF, a12_countCF;
		//case
		int a1_countCa, a2_countCa, a1_homo_countCa, a2_homo_countCa, a12_countCa;
		//case male
		int a1_countCaM, a2_countCaM, a1_homo_countCaM, a2_homo_countCaM, a12_countCaM;
		//case female;
		int a1_countCaF, a2_countCaF, a1_homo_countCaF, a2_homo_countCaF, a12_countCaF;
		//control
		int a1_countCon, a2_countCon, a1_homo_countCon, a2_homo_countCon, a12_countCon;
		//control male
		int a1_countConM, a2_countConM, a1_homo_countConM, a2_homo_countConM, a12_countConM;
		//control female
		int a1_countConF, a2_countConF, a1_homo_countConF, a2_homo_countConF, a12_countConF;

		//groups?
		map<string, vector<int> > gm_allele_counts_o;
		map<string, map<string, int> > gm_geno_counts_o;

		vector<int> m_allele_counts_o, m_allele_counts_om, m_allele_counts_of;
		map<string, int> m_geno_counts_o, m_geno_counts_om, m_geno_counts_of;
		vector<int> m_allele_counts_p, m_allele_counts_pm, m_allele_counts_pf;
	    map<string, int> m_geno_counts_p, m_geno_counts_pm, m_geno_counts_pf;
		vector<int> m_allele_counts_c, m_allele_counts_cm, m_allele_counts_cf;
	    map<string, int> m_geno_counts_c, m_geno_counts_cm, m_geno_counts_cf;
		vector<int> m_allele_counts_ca, m_allele_counts_cam, m_allele_counts_caf;
		map<string, int> m_geno_counts_ca, m_geno_counts_cam, m_geno_counts_caf;
		vector<int> m_allele_counts_con, m_allele_counts_conm, m_allele_counts_conf;
		map<string, int> m_geno_counts_con, m_geno_counts_conm, m_geno_counts_conf;

	public:
		AlleleFrequency(){
			markers = NULL;
			families = NULL;
			samples = NULL;
			threshold = 0.0;
			rank =0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};
		AlleleFrequency(vector<Sample*>* samps){
			samples = samps;
			markers = NULL;
			families = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};
		AlleleFrequency(vector<Sample*>* samps, vector<Family*>* fams){
			samples = samps;
			families = fams;
			markers = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};
		AlleleFrequency(DataSet* ds){
			data_set = ds;
			samples = ds->get_samples();
			families = ds->get_families();
			markers = ds->get_markers();
			marker_map = ds->get_marker_map();
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
		};

		AlleleFrequency(float thresh) : threshold(thresh){
			markers = NULL;
			families = NULL;
			samples = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};
		~AlleleFrequency(){};
//		void process(Connection *, Families*, Markers*);
//		void process(Families*, Markers*);
		void resetDataSet(vector<Sample*>* samps, vector<Family*>* fams){
			samples = samps;
			families = fams;
			markers = NULL;
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			initializeCounts(0);
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
		};
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
		};
		void set_parameters(StepOptions* o){setOptions(*o);};
		void setOptions(StepOptions o){
			options = o;
			if(options.doRandomChild() || options.doAll() || options.doAllChildren() || options.doUnaffSpousesOnly() || options.doUnknownSpouses()){
				flagSamples();
			}
			else{
				options.setFoundersOnly();
			    flagSamples();
			}

		};
		void setOrder(int o){order = o;};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};

		//overall
		int getAonehomo(){return a1_homo_count;};
		int getAtwohomo(){return a2_homo_count;};
		int getHet(){return a12_count;};
		float getAonehomo_exp();
		float getAtwohomo_exp();
		float getHet_exp();
		float getAone_freq();
		float getAtwo_freq();
		int getAone_count(){return a1_count;};
		int getAtwo_count(){return a2_count;};
		int getPop(){return (getAonehomo() + getAtwohomo() + getHet());};
		int getMicroCount(int l){return m_allele_counts_o.at(l);};
		float getMicroFreq(int l);
		int getMicroDenom();
		int getGroupMicroCount(string group, int loc){return gm_allele_counts_o[group][loc];};
		int getGroupAone_count(string group){return ga1_count[group];};
		float getGroupAone_freq(string group){
			if(getGroupPop(group) > 0){
				return ((float)ga1_count[group] / (float)(ga1_count[group] + ga2_count[group]));
			}
			return 0;
		}
		float getGroupAtwo_freq(string group){
			if(getGroupPop(group) > 0){
				return ((float)ga2_count[group] / (float)(ga1_count[group] + ga2_count[group]));
			}
			return 0;
		}
		int getGroupAtwo_count(string group){return ga2_count[group];};
		int getGroupAonehomo(string group){return ga1_homo_count[group];};
		int getGroupHet(string group){return ga12_count[group];};
		int getGroupAtwohomo(string group){return ga2_homo_count[group];};
		int getGroupPop(string group){return (getGroupAonehomo(group) + getGroupHet(group) + getGroupAtwohomo(group));};

		//overall male
		int getAonehomoM(){return a1_homo_countM;};
		int getAtwohomoM(){return a2_homo_countM;};
		int getHetM(){return a12_countM;};
		float getAonehomoM_exp();
		float getAtwohomoM_exp();
		float getHetM_exp();
		float getAoneM_freq();
		float getAtwoM_freq();
		int getAoneM_count(){return a1_countM;};
		int getAtwoM_count(){return a2_countM;};
		int getPopM(){return (getAonehomoM() + getAtwohomoM() + getHetM());};
		int getMicroCountM(int l){return m_allele_counts_om.at(l);};

		//overall female
		int getAonehomoF(){return a1_homo_countF;};
		int getAtwohomoF(){return a2_homo_countF;};
		int getHetF(){return a12_countF;};
		float getAonehomoF_exp();
		float getAtwohomoF_exp();
		float getHetF_exp();
		float getAoneF_freq();
		float getAtwoF_freq();
		int getAoneF_count(){return a1_countF;};
		int getAtwoF_count(){return a2_countF;};
		int getPopF(){return (getAonehomoF() + getAtwohomoF() + getHetF());};
		int getMicroCountF(int l){return m_allele_counts_of.at(l);};

		//parent
		int getAonehomoP(){return a1_homo_countP;};
		int getAtwohomoP(){return a2_homo_countP;};
		int getHetP(){return a12_countP;};
		float getAonehomoP_exp();
		float getAtwohomoP_exp();
		float getHetP_exp();
		float getAoneP_freq();
		float getAtwoP_freq();
		int getAoneP_count(){return a1_countP;};
		int getAtwoP_count(){return a2_countP;};
		int getPopP(){return (getAonehomoP() + getAtwohomoP() + getHetP());};
		int getMicroCountP(int l){return m_allele_counts_p.at(l);};

		//parent MOM
		int getAonehomoPF(){return a1_homo_countPF;};
		int getAtwohomoPF(){return a2_homo_countPF;};
		int getHetPF(){return a12_countPF;};
		float getAonehomoPF_exp();
		float getAtwohomoPF_exp();
		float getHetPF_exp();
		float getAonePF_freq();
		float getAtwoPF_freq();
		int getAonePF_count(){return a1_countPF;};
		int getAtwoPF_count(){return a2_countPF;};
		int getPopPF(){return (getAonehomoPF() + getAtwohomoPF() + getHetPF());};
		int getMicroCountPF(int l){return m_allele_counts_pf.at(l);};

		//parent DAD
		int getAonehomoPM(){return a1_homo_countPM;};
		int getAtwohomoPM(){return a2_homo_countPM;};
		int getHetPM(){return a12_countPM;};
		float getAonehomoPM_exp();
		float getAtwohomoPM_exp();
		float getHetPM_exp();
		float getAonePM_freq();
		float getAtwoPM_freq();
		int getAonePM_count(){return a1_countPM;};
		int getAtwoPM_count(){return a2_countPM;};
		int getPopPM(){return (getAonehomoPM() + getAtwohomoPM() + getHetPM());};
		int getMicroCountPM(int l){return m_allele_counts_pm.at(l);};

		//Child
		int getAonehomoC(){return a1_homo_countC;};
		int getAtwohomoC(){return a2_homo_countC;};
		int getHetC(){return a12_countC;};
		float getAonehomoC_exp();
		float getAtwohomoC_exp();
		float getHetC_exp();
		float getAoneC_freq();
		float getAtwoC_freq();
		int getAoneC_count(){return a1_countC;};
		int getAtwoC_count(){return a2_countC;};
		int getPopC(){return (getAonehomoC() + getAtwohomoC() + getHetC());};

		//Child Male
		int getAonehomoCM(){return a1_homo_countCM;};
		int getAtwohomoCM(){return a2_homo_countCM;};
		int getHetCM(){return a12_countCM;};
		float getAonehomoCM_exp();
		float getAtwohomoCM_exp();
		float getHetCM_exp();
		float getAoneCM_freq();
		float getAtwoCM_freq();
		int getAoneCM_count(){return a1_countCM;};
		int getAtwoCM_count(){return a2_countCM;};
		int getPopCM(){return (getAonehomoCM() + getAtwohomoCM() + getHetCM());};

		//Child Female
		int getAonehomoCF(){return a1_homo_countCF;};
		int getAtwohomoCF(){return a2_homo_countCF;};
		int getHetCF(){return a12_countCF;};
		float getAonehomoCF_exp();
		float getAtwohomoCF_exp();
		float getHetCF_exp();
		float getAoneCF_freq();
		float getAtwoCF_freq();
		int getAoneCF_count(){return a1_countCF;};
		int getAtwoCF_count(){return a2_countCF;};
		int getPopCF(){return (getAonehomoCF() + getAtwohomoCF() + getHetCF());};

		//Case
		int getAonehomoCa(){return a1_homo_countCa;};
		int getAtwohomoCa(){return a2_homo_countCa;};
		int getHetCa(){return a12_countCa;};
		float getAonehomoCa_exp();
		float getAtwohomoCa_exp();
		float getHetCa_exp();
		float getAoneCa_freq();
		float getAtwoCa_freq();
		int getAoneCa_count(){return a1_countCa;};
		int getAtwoCa_count(){return a2_countCa;};
		int getPopCa(){return (getAonehomoCa() + getAtwohomoCa() + getHetCa());};
		int getMicroCountCa(int l){return m_allele_counts_ca.at(l);};

		//Case Male
		int getAonehomoCaM(){return a1_homo_countCaM;};
		int getAtwohomoCaM(){return a2_homo_countCaM;};
		int getHetCaM(){return a12_countCaM;};
		float getAonehomoCaM_exp();
		float getAtwohomoCaM_exp();
		float getHetCaM_exp();
		float getAoneCaM_freq();
		float getAtwoCaM_freq();
		int getAoneCaM_count(){return a1_countCaM;};
		int getAtwoCaM_count(){return a2_countCaM;};
		int getPopCaM(){return (getAonehomoCaM() + getAtwohomoCaM() + getHetCaM());};
		int getMicroCountCaM(int l){return m_allele_counts_cam.at(l);};

		//Case Female
		int getAonehomoCaF(){return a1_homo_countCaF;};
		int getAtwohomoCaF(){return a2_homo_countCaF;};
		int getHetCaF(){return a12_countCaF;};
		float getAonehomoCaF_exp();
		float getAtwohomoCaF_exp();
		float getHetCaF_exp();
		float getAoneCaF_freq();
		float getAtwoCaF_freq();
		int getAoneCaF_count(){return a1_countCaF;};
		int getAtwoCaF_count(){return a2_countCaF;};
		int getPopCaF(){return (getAonehomoCaF() + getAtwohomoCaF() + getHetCaF());};
		int getMicroCountCaF(int l){return m_allele_counts_caf.at(l);};

		//Control
		int getAonehomoCon(){return (a1_homo_countCon);};
		int getAtwohomoCon(){return (a2_homo_countCon);};
		int getHetCon(){return (a12_countCon);};
		float getAonehomoCon_exp();
		float getAtwohomoCon_exp();
		float getHetCon_exp();
		float getAoneCon_freq();
		float getAtwoCon_freq();
		int getPopCon(){return (getAonehomoCon() + getAtwohomoCon() + getHetCon());};
		int getAoneCon_count(){return a1_countCon;};
		int getAtwoCon_count(){return a2_countCon;};
		int getMicroCountCon(int l){return m_allele_counts_con.at(l);};

		//Control Male
		int getAonehomoConM(){return a1_homo_countConM;};
		int getAtwohomoConM(){return a2_homo_countConM;};
		int getHetConM(){return a12_countConM;};
		float getAonehomoConM_exp();
		float getAtwohomoConM_exp();
		float getHetConM_exp();
		float getAoneConM_freq();
		float getAtwoConM_freq();
		int getAoneConM_count(){return a1_countConM;};
		int getAtwoConM_count(){return a2_countConM;};
		int getPopConM(){return (getAonehomoConM() + getAtwohomoConM() + getHetConM());};
		int getMicroCountConM(int l){return m_allele_counts_conm.at(l);};

		//Control Female
		int getAonehomoConF(){return a1_homo_countConF;};
		int getAtwohomoConF(){return a2_homo_countConF;};
		int getHetConF(){return a12_countConF;};
		float getAonehomoConF_exp();
		float getAtwohomoConF_exp();
		float getHetConF_exp();
		float getAoneConF_freq();
		float getAtwoConF_freq();
		int getAoneConF_count(){return a1_countConF;};
		int getAtwoConF_count(){return a2_countConF;};
		int getPopConF(){return (getAonehomoConF() + getAtwohomoConF() + getHetConF());};
		int getMicroCountConF(int l){return m_allele_counts_conf.at(l);};


		void flagSamples();
		void calcOne(Marker*);
		void calcOne(int m){calcOne((*markers).at(m));};
		void calculate(int m){calcOne(m);};
		void calculate(Marker* m){calcOne(m);};
		void calcOneGroups(Marker*);
		void calculateGroups(Marker* m){calcOneGroups(m);};
		void calculateGroups(int m){calcOneGroups((*markers).at(m));};
		void processtest();
		void initializeCounts(int);
		void filterOne(Marker*);
		void setOverwrite(bool v){overwrite = v;};
		Sample* findRandomSample(Sample*, vector<int>&, Marker*);
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
};

#endif
