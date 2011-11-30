#ifndef PROCESSHOMOZYGOUS_H
#define PROCESSHOMOZYGOUS_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
//#include "Markers.h"
//#include "Families.h"
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <General.h>
#include <StepOptions.h>
#include <Homozygous.h>
#include <MethodException.h>
#include <DataSet.h>

using namespace std;
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessHomozygous : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Marker*>* markers;
		vector<Family*>* families;
		vector<int>* marker_map;
		//StepOptions options;
		//Markers* markers;
		//Families* families;
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
		string defaultinsert;
		string projectPath;

		vector<int> homoaffcount;
		vector<int> homounaffcount;
		vector<int> homoallcount;
		vector<int> homominallcount;
		vector<int> homomajallcount;

		vector<Marker*> good_markers;


	public:
		ProcessHomozygous(){
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
		ProcessHomozygous(float thresh) : threshold(thresh){
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
		ProcessHomozygous(string, int, Database*, string);
		~ProcessHomozygous(){};
		//void process(Connection*, Families*, Markers*);
		void PrintSummary();
		void filter();
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
		//void updateFamsMarks(Families* f, Markers* m){
		//	families = f;
		//	markers = m;
		//};
		void setDBOUT(){_DBOUTPUT_ = true;};
		void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		//void process(Families*, Markers*);
		void process(DataSet*);
        void setOrder(int o){order = o;};
		double checkSequence(int, int, int&, int&, double&);
		vector<int> findCommonRoh(Sample*, vector<Marker*>*, vector<int>*);
		void process_wgha();
		vector<Marker*> enabledMarkers();
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
		void perform_homozyg_permutations();
		Sample* findRandomSample(map<Sample*, bool>);
		void run(DataSetObject*);
		void dump2db();
		void create_tables();
		void resize(int);
		void FixOutputName();
};
#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
#endif
