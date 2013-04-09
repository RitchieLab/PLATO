#ifndef EXAMPLEMODULE_H   //change to actual process name
#define EXAMPLEMODULE_H   //change to actual process name

#include <vector>
#include <string>

#include "Process.h"
#include <General.h>
#include <StepOptions.h>
#include <DataSet.h>

//using namespace std;
//using namespace Methods;

class ExampleModule : public ProcessImpl<ExampleModule>{
	private:
		Methods::DataSet* data_set;
		std::vector<Methods::Sample*>* samples;
		std::vector<Methods::Marker*>* markers;
		std::vector<Methods::Family*>* families;
		std::vector<int>* marker_map;
		Methods::StepOptions options;
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
		std::vector<int> homoaffcount;
		std::vector<int> homounaffcount;
		std::vector<int> homoallcount;


	public:
		ExampleModule(){
			data_set = NULL;
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
		};
		ExampleModule(float thresh) : threshold(thresh){
			data_set = NULL;
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
		};
		virtual ~ExampleModule(){
		};
		void PrintSummary();
		void filter();
		void setThreshold(std::string s){
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
		void process(Methods::DataSet*);
        void setOrder(int o){order = o;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		//void run(DataSetObject*){};
		void dump2db(){};
};

#endif
