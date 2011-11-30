#ifndef PROCESSMARKERGENOEFF_H
#define PROCESSMARKERGENOEFF_H

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
#include <DataSet.h>
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <General.h>
#include <StepOptions.h>
#include <MarkerGenoEff.h>

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class ProcessMarkerGenoEff : public Process{
	static string stepname;
	private:
		DataSet* data_set;
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

		vector<int> zeros;
		vector<int> total;
		vector<int> casezeros;
		vector<int> casetotal;
		vector<int> controlzeros;
		vector<int> controltotal;
		map<string, vector<int> > groupzeros;
		map<string, vector<int> > grouptotal;

		vector<Marker*> good_markers;

	public:
		ProcessMarkerGenoEff(){
			data_set = NULL;
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
		ProcessMarkerGenoEff(float thresh) : threshold(thresh){
			data_set = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			_MARKERLIST_ = false;
			_STRATIFY_ = false;
			order = 0;
		};
		~ProcessMarkerGenoEff(){
		zeros.resize(0);
		total.resize(0);
		casezeros.resize(0);
		casetotal.resize(0);
		controlzeros.resize(0);
		controltotal.resize(0);
		};
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
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
#endif
