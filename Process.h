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
		StepOptions options;
		Process(){options.setCovarMissing(opts::_COVAR_MISSING_); options.setTraitMissing(opts::_TRAIT_MISSING_);};
		virtual ~Process();
		//virtual void process(Connection*, Families*, Markers*) = 0;
		//virtual void process(Families*, Markers*) = 0;
		virtual void process(DataSet*) = 0;
		virtual void PrintSummary() = 0;
		virtual void filter() = 0;
		virtual void setThreshold(string) = 0;
		virtual void FilterSummary() = 0;
		virtual void setRank(int) = 0;
		//virtual void updateFamsMarks(Families*, Markers* ) = 0;
		virtual void setDBOUT() = 0;
		virtual void setMarkerList() = 0;
		virtual void setStratify() = 0;
		virtual void setOrder(int) = 0;
		virtual void setOverwrite(bool) = 0;
		virtual bool hasIncExc() = 0;
		StepOptions* getOptions(){return &options;}

};

#endif
