#ifndef STEP_H
#define STEP_H

#include <StepOptions.h>
#include <string>
#include <list>
#include "Process.h"

using namespace std;
using namespace Methods;

class Step{
	private:
		string name;
		string threshold;
		bool multi_thresh;
		int order;
		Process* myprocess;

	public:
		Step(){};
		Step(Process* p) : name(p->getName()), threshold(""), multi_thresh(false), order(0), myprocess(p){}
		Step(string n, string t, bool mt){
			multi_thresh = mt;
			name = n;
			threshold = t;
		};
		virtual ~Step(){
		};

		string getName();
		string getThreshold();
		void setThreshold(string);
		void setOverwrite(bool);
		bool getMultiThresh();
		bool hasIncExc();
		void setMultiThresh(bool v){multi_thresh = v;};
		void setProcess(Process*);
		void setProcess(void*);
		Process* getProcess(){return myprocess;};
		void setName(string s){name = s;};
		void setOrder(int);
		void process(DataSet* ds){myprocess->process(ds);};
		void PrintSummary(){myprocess->PrintSummary();};
		void filter(){myprocess->filter();};

		void FilterSummary(){myprocess->FilterSummary();};
		void setRank(int r){myprocess->setRank(r);};
		void close(){delete myprocess;};
		StepOptions* getOptions(){return myprocess->getOptions();}
};

#endif
