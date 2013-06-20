#ifndef PROCESSMDR_H
#define PROCESSMDR_H


#include "Process.h"

class ProcessMDR : public ProcessImpl<ProcessMDR>{
private:
	static const std::string stepname;

public:
	ProcessMDR(){name="MDR Process";}
	virtual ~ProcessMDR(){}

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
