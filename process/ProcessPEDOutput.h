#ifndef PROCESSPEDOUTPUT_H
#define PROCESSPEDOUTPUT_H

#include "Process.h"

class ProcessPEDOutput : public ProcessImpl<ProcessPEDOutput>{

private:
	static const std::string stepname;

public:
	ProcessPEDOutput(){name="Create PED file";}
	virtual ~ProcessPEDOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
