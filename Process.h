#ifndef PROCESS_H
#define PROCESS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <Globals.h>
//#include "Families.h"
//#include "Markers.h"
#include <Sample.h>
#include <Family.h>
#include <Marker.h>
#include <DataSet.h>
#include <StepOptions.h>

using namespace std;
using namespace Methods;

class Process{

	public:
		Process(){options.setCovarMissing(opts::_COVAR_MISSING_); options.setTraitMissing(opts::_TRAIT_MISSING_);}
		virtual ~Process();

		virtual void process(DataSet*) = 0;
		virtual void PrintSummary() = 0;
		virtual void filter() = 0;
		virtual void setThreshold(string) = 0;
		virtual void FilterSummary() = 0;
		virtual void setRank(int) = 0;
		virtual void setDBOUT() = 0;
		virtual void setMarkerList() = 0;
		virtual void setStratify() = 0;
		virtual void setOrder(int) = 0;
		virtual void setOverwrite(bool) = 0;
		virtual bool hasIncExc() = 0;

		StepOptions* getOptions(){return &options;}
		StepOptions get_options(){return options;}
		void set_options(StepOptions* opts){options = *opts;}
		bool has_results(){return hasresults;}
		vector<string> get_tablename(){return tablename;}
		map<string, vector<string> > get_headers(){return headers;}
		vector<string> get_headers(string table){return headers[table];}
		vector<string> get_primary_table(string table){return primary_table[table];}
		string get_table_nickname(int i){return tablenicknames[i];}
		string get_name(){return name;}
		int get_position(){return position;}
		vector<string> get_filenames(){return filenames;}

	protected:
		StepOptions options;

		string name;
		string batchname;
		int position;
		bool hasresults;
		vector<string> tablename;
		map<string, vector<string> > headers;
		map<string, vector<string> > primary_table;
		vector<string> tablenicknames;
		vector<string> filenames;
};
#endif
