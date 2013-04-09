#ifndef PROCESSALLELEFREQUENCY_H
#define PROCESSALLELEFREQUENCY_H

#include <string>
#include <vector>
#include <map>

#include <Sample.h>
#include <Marker.h>
#include <Family.h>
#include <AlleleFrequency.h>
#include <DataSet.h>
#include <StepOptions.h>

#include "Process.h"

class ProcessAlleleFrequency : public ProcessImpl<ProcessAlleleFrequency>{
	static std::string stepname;
	private:
		Methods::DataSet* data_set;
		std::vector<Methods::Sample*>* samples;
		std::vector<Methods::Marker*>* markers;
		std::vector<Methods::Family*>* families;
		std::vector<int>* marker_map;
		//StepOptions options;
//		Markers* markers;
//		Families* families;
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

		//overall groups
		std::map<std::string, int> ga1_count, ga2_count, ga1_homo_count, ga2_homo_count, ga12_count;

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
		std::map<std::string, std::vector<int> > gm_allele_counts_o;
		std::map<std::string, std::map<std::string, int> > gm_geno_counts_o;

		std::vector<int> m_allele_counts_o, m_allele_counts_om, m_allele_counts_of;
		std::map<std::string, int> m_geno_counts_o, m_geno_counts_om, m_geno_counts_of;
		std::vector<int> m_allele_counts_p, m_allele_counts_pm, m_allele_counts_pf;
	    std::map<std::string, int> m_geno_counts_p, m_geno_counts_pm, m_geno_counts_pf;
		std::vector<int> m_allele_counts_c, m_allele_counts_cm, m_allele_counts_cf;
	    std::map<std::string, int> m_geno_counts_c, m_geno_counts_cm, m_geno_counts_cf;
		std::vector<int> m_allele_counts_ca, m_allele_counts_cam, m_allele_counts_caf;
		std::map<std::string, int> m_geno_counts_ca, m_geno_counts_cam, m_geno_counts_caf;
		std::vector<int> m_allele_counts_con, m_allele_counts_conm, m_allele_counts_conf;
		std::map<std::string, int> m_geno_counts_con, m_geno_counts_conm, m_geno_counts_conf;

	public:
		ProcessAlleleFrequency(){
			data_set = NULL;
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
		ProcessAlleleFrequency(std::vector<Methods::Sample*>* samps){
			data_set = NULL;
			samples = samps;
			markers = NULL;
			families = NULL;
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
		};
		ProcessAlleleFrequency(std::vector<Methods::Sample*>* samps, std::vector<Methods::Family*>* fams){
			data_set = NULL;
			samples = samps;
			families = fams;
			markers = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};
		ProcessAlleleFrequency(Methods::DataSet* ds){
			data_set = ds;
			samples = NULL;
			families = NULL;
			markers = NULL;
			rank = 0;
			order = 0;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
		};
		ProcessAlleleFrequency(float thresh) : threshold(thresh){
			data_set = NULL;
			markers = NULL;
			families = NULL;
			samples = NULL;
			rank = 0;
			order = 0;
		 orig_num_markers = 0;
		 orig_num_families = 0;
		 orig_num_individuals = 0;
		};

		virtual ~ProcessAlleleFrequency(){};
//		void process(Connection *, Families*, Markers*);
//		void process(Families*, Markers*);
		void process(Methods::DataSet*);
		void processtest();
		void PrintSummary();
		void filter();
		void setThreshold(std::string s){
			options.setUp(s);
			//threshold = std::atof(s.c_str());
		};
		void setOptions(Methods::StepOptions o){
			options = o;
		};
		void setOrder(int o){order = o;};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
//		void hw_process(Connection*, Families*, Markers*);
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
//		void updateFamsMarks(Families* f, Markers* m){
//		    families = f;
//		    markers = m;
//		};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};

		void initializeCounts(int);
		void doFilter(Methods::Marker*, Methods::AlleleFrequency*);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};

};

#endif
