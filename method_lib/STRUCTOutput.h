#ifndef STRUCTOUTPUT_H
#define STRUCTOUTPUT_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "Marker.h"
#include "Sample.h"
#include "Family.h"
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"
#include "DataSet.h"

using namespace std;


namespace Methods{
class STRUCTOutput{
	private:
		DataSet* data_set;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		StepOptions options;
		int threshold;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_individuals;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;
		vector<string> filenames;

	public:
		STRUCTOutput(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order =0;
		};
		STRUCTOutput(int thresh) : threshold(thresh){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			samples = NULL;
			marker_map = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		~STRUCTOutput(){};
		void PrintSummary();
		vector<string> get_filenames(){return filenames;}
		void filter();
		void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
		void setThreshold(string s){
			options.setUp(s);
			if(options.getStratFile().length() == 0 && !(options.haveStratification())){
				opts::printLog("STRUCTURE Output step requires the step option of '-strat-file <filenam>' to specify a stratification file.\n");
				throw MethodException("STRUCTURE Output step requires the step option of '-strat-file <filenam>' to specify a stratification file.\n");
			}
			if ((options.haveStratification()))
			{
				cout << "reading stratification file... " << endl;
				options.readStratificationFile(options.getStratFile());
			}
		};
		void setOptions(StepOptions o){
			options = o;
			if(options.getStratFile().length() == 0 && !(options.haveStratification())){
				opts::printLog("STRUCTURE Output step requires the step option of '-strat-file <filename>' to specify a stratification file.\n");
				throw MethodException("STRUCTURE Output step requires the step option of '-strat-file <filename>' to specify a stratification file.\n");
			}
			if ((options.haveStratification()))
			{
				cout << "reading stratification file b ... " << endl;
				options.readStratificationFile(options.getStratFile());
			}
		};
		void calculate(DataSet* ds){
			data_set = ds;
			process(ds->get_samples(), ds->get_families(), ds->get_markers(), ds->get_marker_map());
		};
		void FilterSummary();
		int getOrigNumMarkers(){return orig_num_markers;};
		int getOrigNumFamilies(){return orig_num_families;};
		int getOrigNumIndividuals(){return orig_num_individuals;};
		void setRank(int r){rank = r;};
		int getRank(){return rank;};
        void setOrder(int o){order = o;};
		void setDBOUT(){_DBOUTPUT_ = true;};
		int map_sex(char);
        void setMarkerList(){_MARKERLIST_ = true;};
		void setStratify(){_STRATIFY_ = true;};
		void setOverwrite(bool v){overwrite = v;};
		bool hasIncExc(){return options.doIncExcludedSamples();};
};
};
#endif
