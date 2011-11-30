#ifndef STEP_H
#define STEP_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <StepOptions.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
//#include "MendelianErrors.h"
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
		Step(string n, string t, bool mt){
			multi_thresh = mt;
			name = n;
			threshold = t;
		};
		virtual ~Step(){
		//	if(myprocess){
		//		delete(myprocess);
		//	}
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
		//void process(Connection* con, Families* f, Markers* m) {myprocess->process(con, f, m);};
		//void process(Families* f, Markers* m) {myprocess->process(f, m);};
		void process(DataSet* ds){myprocess->process(ds);};
		//Families* f, Markers* m) {myprocess->process(f, m);};
		void PrintSummary(){myprocess->PrintSummary();};
		void filter(){myprocess->filter();};

		void FilterSummary(){myprocess->FilterSummary();};
		void setRank(int r){myprocess->setRank(r);};
		void close(){delete myprocess;};
		StepOptions* getOptions(){return myprocess->getOptions();}
		//void updateFamsMarks(Families* f, Markers* m){myprocess->updateFamsMarks(f, m);};
};

#endif
