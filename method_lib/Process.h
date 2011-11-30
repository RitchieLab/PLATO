#ifndef PROCESS_H
#define PROCESS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "Globals.h"
#include "Sample.h"
#include "Family.h"
#include "Marker.h"
using namespace std;
class Process{
	public:
		Process(){};
		virtual ~Process();
		virtual void process(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*) = 0;
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
};

#endif
