#ifndef PROCESSHWEQUILIBRIUM_H
#define PROCESSHWEQUILIBRIUM_H

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
#include <Marker.h>
#include <Family.h>
#include <Sample.h>
#include <Globals.h>
#include "Process.h"
#include <Options.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <HWEquilibrium.h>
#include <MethodException.h>
#include <AlleleFrequency.h>

using namespace std;
using namespace Methods;

class ProcessHWEquilibrium : public Process{
	static string stepname;
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
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
		bool useoverall;
		bool overwrite;
		int order;

		float hw_O;
		float hw_OM;
		float hw_OF;
		float hw_P;
		float hw_PM;
		float hw_PF;
		float hw_C;
		float hw_CM;
		float hw_CF;
		float hw_Ca;
		float hw_CaF;
		float hw_CaM;
		float hw_Con;
		float hw_ConM;
		float hw_ConF;

		vector<float> hwO;
		vector<float> hwP;
		vector<float> hwPM;
		vector<float> hwPD;
		vector<float> hwC;
		vector<float> hwCM;
		vector<float> hwCF;
		vector<float> hwCa;
		vector<float> hwCaF;
		vector<float> hwCaM;
		vector<float> hwCon;
		vector<float> hwConF;
		vector<float> hwConM;


	public:
		ProcessHWEquilibrium(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			//af = new AlleleFrequency();
			rank = 0;
			_MARKERLIST_ = false;
			_DBOUTPUT_ = false;
			_STRATIFY_ = false;
			order = 0;
			orig_num_markers = 0;
		};
		ProcessHWEquilibrium(float thresh) : threshold(thresh){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			//af = new AlleleFrequency();
			rank = 0;
			_MARKERLIST_ = false;
			_DBOUTPUT_ = false;
			_STRATIFY_ = false;
			order = 0;
			orig_num_markers = 0;
		};
		~ProcessHWEquilibrium(){
			//if(af){
			//delete(af);
		//	}
		};
//		void process(Connection*, Families*, Markers*);
//		void process(Families*, Markers*);
		void PrintSummary();
		void filter();
		void process(DataSet*);
		void setThreshold(string s){
			options.setUp(s);
			//	threshold = std::atof(s.c_str());
		};
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

		void resize(int i);
		void doFilter(Marker*, HWEquilibrium*);
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};

};

#endif
