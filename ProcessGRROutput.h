#ifndef PROCESSGRROUTPUT_H
#define PROCESSGRROUTPUT_H

#include "Process.h"

class ProcessGRROutput : public ProcessImpl<ProcessGRROutput>{
private:
	static const std::string stepname;

public:
	ProcessGRROutput(){name="Create GRR input file";}
	virtual ~ProcessGRROutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
