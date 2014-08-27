#ifndef PROCESSMDRPDT_H
#define PROCESSMDRPDT_H

#include "Process.h"

class ProcessMDRPDT : public ProcessImpl<ProcessMDRPDT>{

private:
	static const std::string stepname;

public:
	ProcessMDRPDT(){name="MDRPDT Process";}
	virtual ~ProcessMDRPDT(){}

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);

private:
	void doFilter(Methods::Marker*, double);
};
#endif
