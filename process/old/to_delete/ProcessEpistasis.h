#ifndef PROCESSEPISTASIS_H
#define PROCESSEPISTASIS_H

#include "Process.h"

class ProcessEpistasis : public ProcessImpl<ProcessEpistasis>{

private:
	static const std::string stepname;

public:
	ProcessEpistasis(){name="Epistasis";}
	virtual ~ProcessEpistasis(){};

	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
