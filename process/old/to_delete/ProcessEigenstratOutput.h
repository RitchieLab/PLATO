#ifndef PROCESSEIGENSTRATOUTPUT_H
#define PROCESSEIGENSTRATOUTPUT_H

#include "Process.h"

class ProcessEigenstratOutput : public ProcessImpl<ProcessEigenstratOutput>{

private:
	static const std::string stepname;

public:
	ProcessEigenstratOutput(){name="Create Eigenstrat file";}
	virtual ~ProcessEigenstratOutput(){};

protected:
	virtual void process(Methods::DataSet*);
	void PrintSummary();
};
#endif
