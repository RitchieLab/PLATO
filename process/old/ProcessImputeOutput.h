#ifndef PROCESSIMPUTEOUTPUT_H
#define PROCESSIMPUTEOUTPUT_H

#include "Process.h"


using namespace std;
using namespace Methods;

class ProcessImputeOutput : public ProcessImpl<ProcessImputeOutput>{

private:
	static const string stepname;

public:
	ProcessImputeOutput(){name="IBS";}
	virtual ~ProcessImputeOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
