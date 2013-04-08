#ifndef PERCENTBYFAMILY_H
#define PERCENTBYFAMILY_H

#include <vector>
#include <string>
#include <map>

#include "Process.h"
#include <General.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <Marker.h>

//typedef vector<int> PERCENT;

class PercentByFamily : public Process{
	static std::string stepname;
	private:
		Methods::DataSet* data_set;
		std::vector<Methods::Sample*>* samples;
		std::vector<Methods::Marker*>* markers;
		std::vector<Methods::Family*>* families;
		std::vector<int>* marker_map;
		StepOptions options;
		//Markers* markers;
		//Families* families;
		float threshold;
		//PERCENT per_cutoff;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_samples;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;

		std::vector<int> fzeros;
		std::vector<int> ftotal;
		std::vector<int> mzeros;
		std::vector<int> mtotal;
		std::vector< std::map<std::string, int> > enzyme_zeros;
		std::vector< std::map<std::string, int> > enzyme_total;

	public:
		PercentByFamily(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			rank = 0;
			order = 0;
		};
		PercentByFamily(float thresh) : threshold(thresh){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_samples = 0;
			rank = 0;
			order =0;
		};
		virtual ~PercentByFamily(){};


//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(Methods::DataSet*);
		void setThreshold(std::string s){threshold = std::atof(s.c_str());};
		//PERCENT* getPerCutoff(){return &per_cutoff;};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_samples;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
//		void updateFamsMarks(Families* f, Markers* m){
//			families = f;
//			markers = m;
//		};
        void setDBOUT(){_DBOUTPUT_ = true;};
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
        void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		//void run(DataSetObject*){};
		void dump2db(){};
};

#endif
