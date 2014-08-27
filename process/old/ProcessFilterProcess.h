#ifndef PROCESSFILTERPROCESS_H
#define PROCESSFILTERPROCESS_H


#include "Process.h"

using namespace std;
using namespace Methods;

class ProcessFilterProcess : public ProcessImpl<ProcessFilterProcess>{
private:
	static const string stepname;

public:
	ProcessFilterProcess(){name="Filter Process";}
	virtual ~ProcessFilterProcess(){};

protected:
	virtual	void PrintSummary();
	virtual void process(Methods::DataSet*);
};

#endif
