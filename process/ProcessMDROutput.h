#ifndef PROCESSMDROUTPUT_H
#define PROCESSMDROUTPUT_H

#include "Process.h"

class ProcessMDROutput : public ProcessImpl<ProcessMDROutput>{

private:
	static const std::string stepname;

public:
	ProcessMDROutput(){name="Create MDR input file";}
	virtual ~ProcessMDROutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
