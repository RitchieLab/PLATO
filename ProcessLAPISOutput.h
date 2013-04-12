#ifndef PROCESSLAPISOUTPUT_H
#define PROCESSLAPISOUTPUT_H

#include "Process.h"

class ProcessLAPISOutput : public ProcessImpl<ProcessLAPISOutput>{

private:
	static const std::string stepname;

public:
	ProcessLAPISOutput(){name="Create LAPIS input file";}
	virtual ~ProcessLAPISOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
